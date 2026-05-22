#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# PROGRAM NAME: iedb_scraper.py
# PURPOSE:
#   Read a file with IEDB links, download each page with httpx,
#   and extract:
#     1) epitope
#     2) ALL antigen + UniProt + organism pairs from referenceEpitopeString
#     3) MHC ligand summary (positive and negative alleles)
#     4) T cell assay summary (total response + assay type columns)
#   Then save:
#     - a WIDE CSV  -> easy to open in Excel
#     - a LONG CSV  -> easier for analysis in pandas / R
# -----------------------------------------------------------------------------

# ------------------------------- Libraries -----------------------------------

import argparse  # Read options from the terminal, like -i and -o
import json  # Convert JSON text into Python dictionaries
import logging  # Save progress messages and errors
import re  # Search text with regular expressions
import time  # Lets Python pause between requests with time.sleep()
import random  # Lets Python add a small random wait so requests are less mechanical
from collections import OrderedDict  # Keep dictionary keys in a nice order
from pathlib import Path  # Work with file paths more safely

import httpx  # Download HTML pages from the web
import pandas as pd  # Build and save tables
from alive_progress import alive_bar  # Show a progress bar while links are processed
from bs4 import BeautifulSoup  # Parse HTML and inspect <script> tags
import csv

# ----------------------------- Helper functions ------------------------------


def parse_args():
    """
    Read arguments from the terminal.

    Example:
        python script.py -i links.txt -o output.csv
        python script.py -i links.txt -o output.csv --output_long_csv output_long.csv
    """
    parser = argparse.ArgumentParser(
        description="Parse IEDB HTML pages and generate wide and long CSV tables."
    )

    # Input file: one IEDB link per line
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Text file with one IEDB link per line.",
    )

    # Output CSV file for the wide table
    parser.add_argument(
        "-o",
        "--output_csv",
        required=True,
        help="Output CSV file name for the wide table.",
    )

    # Optional output CSV file for the long table
    parser.add_argument(
        "-o2",
        "--output_long_csv",
        default=None,
        help=(
            "Optional output CSV file name for the long table. "
            'If omitted, the script will create one automatically, for example "results_long.csv".'
        ),
    )

    # Optional organism override, kept from the old script
    parser.add_argument(
        "-org",
        "--organism",
        help="Optional organism override for the WIDE table.",
    )

    # Optional log file
    parser.add_argument(
        "-l",
        "--log",
        default="log.txt",
        help='Log file name (optional). Default is "log.txt".',
    )

    # Optional normal sleep time between links.
    # Example:
    #   -s 2
    # means:
    #   after finishing one link, wait 2 seconds before starting the next link.
    #
    # Use 0 if you want to go as fast as possible between successful links.
    parser.add_argument(
        "-s",
        "--sleep_seconds",
        type=float,
        default=0.0,
        help=(
            "Number of seconds to wait BETWEEN normal links. "
            "Use 0 to avoid normal waiting. Default: 0."
        ),
    )

    # Optional base delay for 429 retry waits.
    # Example:
    #   -r 20
    # means:
    #   if IEDB returns 429 and does NOT send Retry-After, start by waiting ~20 seconds.
    #
    # This is separate from --sleep_seconds on purpose.
    # --sleep_seconds controls normal waiting between different links.
    # --retry_base_delay controls emergency waiting when the server says 429.
    parser.add_argument(
        "-r",
        "--retry_base_delay",
        type=float,
        default=20.0,
        help=(
            "Base number of seconds to wait before retrying after a 429 error "
            "when the server does not provide Retry-After. Default: 20."
        ),
    )

    # the script can try the same link up to 5 times when it gets a retryable error like 429
    parser.add_argument(
        "-m",
        "--max_retries",
        type=int,
        default=5,
        help="Maximum number of download attempts per link. Default: 5.",
    )

    return parser.parse_args()


def get_backoff_wait_seconds(retry_base_delay, attempt):
    """
    Calculate how long to wait before retrying a failed download.

    Child version:
        The first time the robot fails, it waits a little.
        If it fails again, it waits longer.
        If it fails again, it waits even longer.

    This is called exponential backoff.

    Example with retry_base_delay = 20:
        attempt 1 -> about 20 seconds
        attempt 2 -> about 40 seconds
        attempt 3 -> about 80 seconds

    We also add 0 to 3 random seconds so the requests do not happen
    at perfectly robotic times.
    """
    return retry_base_delay * (2 ** (attempt - 1)) + random.uniform(0, 3)


def get_long_output_path(output_csv, output_long_csv=None):
    """
    Decide the long output file name.

    If the user gave --output_long_csv, use that.
    Otherwise:
        results.csv -> results_long.csv
    """
    if output_long_csv:
        return output_long_csv

    output_path = Path(output_csv)
    long_name = f"{output_path.stem}_long{output_path.suffix}"
    return str(output_path.with_name(long_name))


def get_failed_output_path(output_csv):
    """
    Decide where to save the failed-links CSV.

    Important idea:
    - We want the failed CSV to live in the SAME DIRECTORY as the main output CSV.
    - That way, all output files from one run stay together.

    Example:
        output_csv = /path/to/results.csv
        failed CSV = /path/to/results_failed_links.csv

    Child version:
        If the main output goes into one folder, the failed-links list
        should go into the same folder, like putting all school papers
        into the same backpack.
    """
    output_path = Path(output_csv)
    failed_name = f"{output_path.stem}_failed_links{output_path.suffix}"
    return str(output_path.with_name(failed_name))


def get_iedb_link_list(file_path):
    """
    Open the input file and return a list of links.
    """
    links_list = []

    with open(file_path, mode="r", encoding="utf-8") as fh:
        for link in fh:
            link = link.strip()
            if link:
                links_list.append(link)

    return links_list


def normalize_iedb_url(url):
    """
    Normalize IEDB URLs so they use https://

    Why do this?
    - Old link lists sometimes still use http://
    - We want one consistent format in the final tables
    """
    url = url.strip()

    if url.startswith("http://"):
        url = "https://" + url[len("http://") :]

    return url


def build_httpx_client():
    """
    Build one reusable HTTPX client.

    Why use one client?
    - It can reuse connections.
    - It keeps all headers in one place.
    - In your testing, httpx worked where requests did not.
    """
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/124.0 Safari/537.36"
        ),
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9",
        "Referer": "https://www.iedb.org/",
    }

    return httpx.Client(headers=headers, follow_redirects=True, timeout=30.0)


def extract_html(url, client, retry_base_delay=20.0, max_retries=5):
    """
    Download the HTML of one webpage and return it as text.

    This function is the DOWNLOADER.

    Simple explanation:
        The script is like a little robot knocking on IEDB's door.
        Usually, the door opens and the robot gets the page.
        But sometimes IEDB says:
            "429 Too Many Requests. You are knocking too fast."
        There can be other temporary problems:
            - server disconnected without sending a response
            - timeout
            - temporary server error
        When that happens, this function waits and tries the SAME link again.

    Parameters
    ----------
    url : str
        The IEDB page URL.

    client : httpx.Client
        The reusable HTTPX client that downloads pages.

    retry_base_delay : float
        The base number of seconds to wait after a 429 error IF the server
        does not provide a Retry-After header.

        Example:
            retry_base_delay = 20
            first 429 retry  -> wait about 20 seconds
            second 429 retry -> wait about 40 seconds
            third 429 retry  -> wait about 80 seconds

    max_retries : int
        How many times to try the same URL before giving up.
    """
    # Make sure the URL uses https:// instead of http://.
    url = normalize_iedb_url(url)

    # These HTTP status codes often mean:
    #   "The problem may be temporary; waiting and retrying can help."
    # We do NOT retry all errors. For example, 404 usually means the page
    # does not exist, so retrying many times would not help.
    retryable_status_codes = {429, 500, 502, 503, 504}
    # 429 means: Too Many Requests.
    # The server is asking us to slow down.

    # Try the same URL multiple times.
    # attempt will be 1, 2, 3, ..., max_retries.
    for attempt in range(1, max_retries + 1):
        try:
            # Ask IEDB for the page.
            response = client.get(url)

            # If the server returned an error like 429, 404, or 500,
            # this line raises an HTTPStatusError.
            # this line raises an exception.
            response.raise_for_status()

            # If we reached this line, the download worked.
            # We return the HTML text to the parser.
            return response.text

        except httpx.HTTPStatusError as e:
            # HTTPStatusError means the server answered, but with
            # an error status code such as 429, 403, 404, or 500.
            status_code = e.response.status_code

            # Only retry errors that are likely to be temporary.
            if status_code in retryable_status_codes and attempt < max_retries:
                # Some servers send a Retry-After header.
                # Example:
                #   Retry-After: 60
                # means:
                #   wait 60 seconds before trying again.
                retry_after = e.response.headers.get("Retry-After")

                if retry_after and retry_after.isdigit():
                    # The server gave us an exact waiting time.
                    # The server's instruction wins over our local setting.
                    wait_seconds = float(retry_after)
                    wait_reason = "server Retry-After header"
                else:
                    # The server did NOT tell us exactly how long to wait.
                    # So we use our local retry_base_delay argument.
                    wait_seconds = get_backoff_wait_seconds(retry_base_delay, attempt)
                    wait_reason = "local retry_base_delay"

                logging.warning(
                    f"HTTP {status_code} for {url}. "
                    f"Waiting {wait_seconds:.1f} seconds before retry "
                    f"{attempt + 1}/{max_retries} "
                    f"using {wait_reason}."
                )

                # Make Python sleep before trying the SAME URL again.
                time.sleep(wait_seconds)
                continue

            # If the error is not retryable, or if we already used
            # used all retries, send the error upward to main().
            raise

        except httpx.RequestError as e:
            # RequestError means the request failed before we got a normal
            # HTTP response. This includes cases like:
            #   - Server disconnected without sending a response
            #   - timeout
            #   - network connection problem
            #
            # These problems are often temporary, so retrying can help.
            if attempt < max_retries:
                wait_seconds = get_backoff_wait_seconds(retry_base_delay, attempt)

                logging.warning(
                    f"Network/protocol error for {url}: {e}. "
                    f"Waiting {wait_seconds:.1f} seconds before retry "
                    f"{attempt + 1}/{max_retries}."
                )

                # Make Python sleep before trying the SAME URL again.
                time.sleep(wait_seconds)
                continue

            # If all retries failed, send the error upward to main().
            raise

    # Safety net. Normally Python reaches this only if something unusual happens.
    raise RuntimeError(f"Could not download {url} after {max_retries} attempts.")


def clean_whitespace(text):
    """
    Replace repeated spaces/newlines/tabs with single spaces.
    """
    return " ".join(text.split())


def get_reference_epitope_string(script_text):
    """
    Find the JavaScript variable called 'refernceEpitopeData'
    and return the value of:
        data -> referenceEpitopeString

    Important:
    The IEDB page uses the misspelled name 'refernceEpitopeData',
    so our code must use that exact spelling.
    """
    pattern = r"var\s+refernceEpitopeData\s*=\s*(\{.*?\});"
    match = re.search(pattern, script_text, flags=re.DOTALL)

    if not match:
        raise ValueError("Could not find refernceEpitopeData in this script block.")

    json_text = match.group(1)
    ref_obj = json.loads(json_text)

    try:
        reference_string = ref_obj["data"]["referenceEpitopeString"]
    except (KeyError, TypeError):
        raise ValueError("referenceEpitopeString was not found in refernceEpitopeData.")

    return reference_string


def extract_epitope(reference_string):
    """
    Extract the epitope from the beginning of the sentence.

    Example:
        'AAISDYDYY is a linear peptidic epitope ...'
    becomes:
        'AAISDYDYY'
    """
    epitope_match = re.search(r"^(.*?)\s+is\s+(?:a|an)\b", reference_string)

    if not epitope_match:
        raise ValueError("Could not extract the epitope from referenceEpitopeString.")

    return clean_whitespace(epitope_match.group(1))


def extract_antigen_organism_pairs(reference_string):
    """
    Extract ALL antigen-organism pairs from the reference sentence.

    This function looks inside a sentence like:

        KLWAQCVQL is a linear peptidic epitope ...
        studied as part of
        Replicase polyprotein 1ab (UniProt:P0C6X7) from SARS-CoV1,
        Replicase polyprotein 1ab (UniProt:P0DTD1) from SARS-CoV2,
        Replicase polyprotein 1a (UniProt:P0DTC1) from SARS-CoV2
        and Replicase polyprotein 1a (UniProt:P0C6U8) from SARS-CoV1.

    and extracts this:

        Pair 1:
            antigen  = Replicase polyprotein 1ab
            uniprot  = P0C6X7
            organism = SARS-CoV1

        Pair 2:
            antigen  = Replicase polyprotein 1ab
            uniprot  = P0DTD1
            organism = SARS-CoV2

        Pair 3:
            antigen  = Replicase polyprotein 1a
            uniprot  = P0DTC1
            organism = SARS-CoV2

        Pair 4:
            antigen  = Replicase polyprotein 1a
            uniprot  = P0C6U8
            organism = SARS-CoV1

    Important:
    - It works with pairs separated by commas.
    - It works with the last pair separated by 'and'.
    - It can handle missing UniProt codes.
    - It does not care whether there are 1, 2, 3, 4, or more pairs.
    """

    # ------------------------------------------------------------
    # Step 1: isolate the useful part of the sentence
    # ------------------------------------------------------------
    #
    # Full reference string example:
    #
    #   "KLWAQCVQL is ... studied as part of Protein A from Virus 1,
    #    Protein B from Virus 2 and Protein C from Virus 3.
    #    This epitope has been studied ..."
    #
    # We only want the part after:
    #
    #   studied as part of
    #
    # and before:
    #
    #   . This epitope
    #
    # So we extract only:
    #
    #   Protein A from Virus 1, Protein B from Virus 2 and Protein C from Virus 3
    # ------------------------------------------------------------

    section_match = re.search(
        # Stop at the start of the next epitope sentence.
        #
        # Why this matters:
        # Some pages say:
        #   ". This epitope has been studied ..."
        # Other pages say:
        #   ". The epitope is an analog of ..."
        #
        # If we do not stop before both options, the organism can accidentally
        # include extra text like:
        #   "SARS-CoV2. The epitope is an analog of ..."
        # r"studied as part of\s+(.*?)(?=\.\s+This epitope\b|$)",
        r"studied as part of\s+(.*?)(?=\.\s+(?:This|The)\s+epitope\b|$)",
        reference_string,
        flags=re.DOTALL,
    )

    # If we cannot find the "studied as part of" section,
    # we stop and tell the user/script that parsing failed.
    if not section_match:
        raise ValueError("Could not isolate the antigen-organism section.")

    # Clean extra spaces and line breaks.
    pairs_section = clean_whitespace(section_match.group(1))

    # Sometimes the extracted section may end with a final period.
    # We remove it so the final organism does not become:
    #   "SARS-CoV1."
    # but instead:
    #   "SARS-CoV1"
    pairs_section = pairs_section.rstrip(".")

    # ------------------------------------------------------------
    # Step 2: define a regex that can find MANY pairs
    # ------------------------------------------------------------
    #
    # This regex understands that pairs may start:
    #
    #   1. at the beginning of the section
    #   2. after a comma
    #   3. after the word "and"
    #
    # Examples it can read:
    #
    #   Protein A (UniProt:AAA) from Virus 1,
    #   Protein B (UniProt:BBB) from Virus 2
    #   and Protein C (UniProt:CCC) from Virus 3
    #
    # It also works if UniProt is missing:
    #
    #   Other SARS-CoV1 protein from SARS-CoV1
    # ------------------------------------------------------------

    pair_pattern = re.compile(
        r"(?:^|,\s+|\s+and\s+)"  # pair starts at beginning, after comma, or after 'and'
        r"(?P<antigen_name>.+?)"  # capture antigen/protein name
        r"(?:\s+\((?P<uniprot>UniProt:[A-Za-z0-9_-]+)\))?"  # optionally capture UniProt code
        r"\s+from\s+"  # match the word 'from'
        r"(?P<organism>.*?)"  # capture organism name
        r"(?=(?:,\s+|\s+and\s+).+?(?:\s+\(UniProt:[A-Za-z0-9_-]+\))?\s+from\s+|$)",  # stop before next pair or end
        flags=re.DOTALL,
    )

    # This list will store all extracted pairs.
    pairs = []

    # ------------------------------------------------------------
    # Step 3: find all matches
    # ------------------------------------------------------------
    #
    # finditer() means:
    #   "Find every match, not just the first one."
    #
    # This is why the function can extract 3, 4, 5, or more pairs.
    # ------------------------------------------------------------

    for match in pair_pattern.finditer(pairs_section):
        # Extract and clean the antigen name.
        antigen_name = clean_whitespace(match.group("antigen_name"))

        # Extract and clean the organism name.
        organism_name = clean_whitespace(match.group("organism"))

        # Extract the UniProt part.
        # It may be something like:
        #   UniProt:P0C6X7
        #
        # Or it may be None if no UniProt code was present.
        raw_uniprot = match.group("uniprot")

        if raw_uniprot:
            # Convert:
            #   UniProt:P0C6X7
            # into:
            #   P0C6X7
            uniprot_code = raw_uniprot.replace("UniProt:", "")

            # A human-readable combined version.
            antigen_full = f"{antigen_name} (UniProt:{uniprot_code})"
        else:
            # If there is no UniProt code, store "-".
            uniprot_code = "-"
            antigen_full = antigen_name

        # Save this one antigen-UniProt-organism set.
        pairs.append(
            {
                "antigen": antigen_name,
                "uniprot": uniprot_code,
                "organism": organism_name,
                "antigen_full": antigen_full,
            }
        )

    # ------------------------------------------------------------
    # Step 4: fail clearly if nothing was extracted
    # ------------------------------------------------------------
    #
    # This is useful for debugging.
    # If a future IEDB sentence has a strange structure,
    # the log will show the exact section that failed.
    # ------------------------------------------------------------

    if not pairs:
        logging.error(f"Could not extract pairs from pairs_section: {pairs_section}")
        raise ValueError("Could not extract any antigen-organism pair.")

    return pairs


def get_compiled_data(script_text):
    """
    Find the JavaScript variable called 'compiledData' and convert it to Python.
    """
    pattern = r"var\s+compiledData\s*=\s*(\{.*?\});"
    match = re.search(pattern, script_text, flags=re.DOTALL)

    if not match:
        raise ValueError("Could not find compiledData in this script block.")

    try:
        return json.loads(match.group(1))
    except json.JSONDecodeError:
        raise ValueError("Failed to parse compiledData as JSON.")


def get_compiled_section_by_type(compiled_data_json, target_type):
    """
    Look inside compiledData['data'] and return the section whose 'type'
    matches target_type.

    This is safer than assuming:
    - first section = MHC
    - last section = T cell

    because we are using the section names directly.
    """
    try:
        sections = compiled_data_json.get("data", [])
    except AttributeError:
        return []

    for section in sections:
        if section.get("type") == target_type:
            return section.get("data", [])

    return []


def allele_classification(mhc_list):
    """
    Split MHC molecules into positive and negative lists.

    Rule:
    - if positive_count > 0 -> positive list
    - else -> negative list
    """
    mhc_positives = []
    mhc_negatives = []

    for mhc_dict in mhc_list:
        if "mhc_molecule" in mhc_dict and "positive_count" in mhc_dict:
            try:
                mhc_molecule = mhc_dict["mhc_molecule"]
                mhc_positive_count = int(mhc_dict["positive_count"])

                if mhc_positive_count > 0:
                    mhc_positives.append(mhc_molecule)
                else:
                    mhc_negatives.append(mhc_molecule)
            except ValueError:
                continue

    return mhc_positives, mhc_negatives


def build_tcell_assay_info_list(tcell_list):
    """
    Convert the raw list of T cell assay dictionaries into a list of tuples.

    Example:
        [("qualitative binding", "1/1")]
    """
    assay_info_list = []

    for assay_dict in tcell_list:
        assay_type = assay_dict.get("assay_type", "Unknown assay type")
        positive_count = assay_dict.get("positive_count", "0")
        total_count = assay_dict.get("total_count", "0")
        assay_result = f"{positive_count}/{total_count}"
        assay_info_list.append((assay_type, assay_result))

    return assay_info_list


def calculate_total_t_response(assay_list):
    """
    Calculate the total T cell response.
    i.e.:
    Calculates if the total T cell response is positive (1) or negative (0).

    Rule used in the old script:
    - if any assay has a numerator > 0, return 1: response is positive
    - otherwise return 0: response is negative
    - if there are no assays, return '-'
    """
    if not assay_list:
        return "-"

    t_assays_boolean_list = [int(val.split("/")[0]) > 0 for _, val in assay_list]
    return 1 if any(t_assays_boolean_list) else 0


def extract_all_data_from_html(html, source_link, failed_csv):
    """
    From one HTML page, find the right <script> block and extract:
    - epitope
    - all antigen-organism-UniProt pairs
    - MHC positive/negative allele lists
    - T cell assay summary
    - source link
    """
    soup = BeautifulSoup(html, "html.parser")
    scripts = soup.find_all("script", type="text/javascript")

    # This flag remembers whether we found a script block that looked like
    # the IEDB detail-data block we need.
    #
    # Why?
    # - If we never find the expected script block, that is also a parsing failure.
    # - We want to save that URL into the failed-links CSV.
    # - But if we did find the block and parsing failed inside it, that failure
    #   is already saved in the except blocks below, so we should not save it twice.
    found_candidate_script = False

    for script in scripts:
        script_text = script.string if script.string else script.get_text()

        if not script_text:
            continue

        # We need BOTH variables in the same script block.
        if (
            "refernceEpitopeData" not in script_text
            or "compiledData" not in script_text
        ):
            continue

        # If we reached this line, we found the script block that probably
        # contains the IEDB detail data.
        found_candidate_script = True

        # Start with a placeholder.
        # If parsing fails before we extract the reference string,
        # we can still write something safe to the failed CSV.
        reference_string = "-"

        try:
            reference_string = get_reference_epitope_string(script_text)
            epitope = extract_epitope(reference_string)
            pairs = extract_antigen_organism_pairs(reference_string)

            compiled_data_json = get_compiled_data(script_text)

            mhc_list = get_compiled_section_by_type(compiled_data_json, "mhc_ligand")
            positive_list, negative_list = allele_classification(mhc_list)

            tcell_list = get_compiled_section_by_type(compiled_data_json, "tcell")
            assay_info_list = build_tcell_assay_info_list(tcell_list)
            total_t_response = calculate_total_t_response(assay_info_list)

            assay_entries = OrderedDict()
            for assay_key, assay_val in assay_info_list:
                # Keep the first occurrence of each assay type column
                if assay_key not in assay_entries:
                    assay_entries[assay_key] = assay_val

            return {
                "Epitope": epitope,
                "pairs": pairs,
                "Positive MHC alleles": (
                    ",".join(positive_list) if positive_list else "-"
                ),
                "Negative MHC alleles": (
                    ",".join(negative_list) if negative_list else "-"
                ),
                "Total response T cell assay(s)": total_t_response,
                "assay_entries": assay_entries,
                "Source": source_link,
            }

        except ValueError as e:
            logging.exception(f"ValueError while parsing {source_link}: {e}")

            # Save parsing failures to the failed-links CSV.
            #
            # Example parsing failure:
            #   The HTML was downloaded, but the script could not extract
            #   antigen-organism pairs from the reference sentence.
            save_failed_link(
                failed_csv=failed_csv,
                source_link=source_link,
                error_message=str(e),
                failure_type="parsing",
                reference_string=reference_string,
            )
        except Exception as e:
            logging.exception(f"Unexpected error while parsing {source_link}: {e}")

            # Save unexpected parser failures too.
            save_failed_link(
                failed_csv=failed_csv,
                source_link=source_link,
                error_message=str(e),
                failure_type="parsing_unexpected",
                reference_string=reference_string,
            )

    # If we never found the expected IEDB script block, save this as a
    # parsing failure too.
    #
    # Child version:
    #   The downloader brought us a book, but the book did not contain
    #   the page/section that our parser expected to read.
    if not found_candidate_script:
        save_failed_link(
            failed_csv=failed_csv,
            source_link=source_link,
            error_message=(
                "Could not find a script block containing both "
                "refernceEpitopeData and compiledData."
            ),
            failure_type="parsing_no_expected_script",
            reference_string="-",
        )

    return {
        "Epitope": "-",
        "pairs": [],
        "Positive MHC alleles": "-",
        "Negative MHC alleles": "-",
        "Total response T cell assay(s)": "-",
        "assay_entries": OrderedDict(),
        "Source": source_link,
    }


def collect_all_assay_column_names(scraped_rows):
    """
    Collect every assay type seen across all rows.

    Why?
    - The old script created one output column per assay type.
    - Different epitopes may have different assay type names.
    """
    assay_names = []

    for row in scraped_rows:
        for assay_name in row["assay_entries"].keys():
            if assay_name not in assay_names:
                assay_names.append(assay_name)

    return assay_names


def flatten_rows(scraped_rows, assay_column_names, organism_override=None):
    """
    Convert scraped rows into a WIDE table.

    Why wide?
    - Easy to open in Excel.
    - Good when most epitopes have only one antigen-organism pair.

    Example columns:
        Epitope,
        Antigen_1, UniProt_1, Organism_1,
        Antigen_2, UniProt_2, Organism_2,
        Positive MHC alleles,
        Negative MHC alleles,
        Total response T cell assay(s),
        qualitative binding,
        Source
    """
    max_pairs = max((len(row["pairs"]) for row in scraped_rows), default=0)

    if max_pairs == 0:
        max_pairs = 1

    flat_rows = []

    for row in scraped_rows:
        flat_row = OrderedDict()
        flat_row["Epitope"] = row["Epitope"]

        # Add numbered antigen / UniProt / organism columns.
        for i in range(1, max_pairs + 1):
            if i <= len(row["pairs"]):
                pair = row["pairs"][i - 1]
                flat_row[f"Antigen_{i}"] = pair["antigen"]
                flat_row[f"UniProt_{i}"] = pair["uniprot"]
                flat_row[f"Organism_{i}"] = organism_override or pair["organism"]
            else:
                flat_row[f"Antigen_{i}"] = "-"
                flat_row[f"UniProt_{i}"] = "-"
                flat_row[f"Organism_{i}"] = "-"

        # Add the old assay-related summary columns.
        flat_row["Positive MHC alleles"] = row["Positive MHC alleles"]
        flat_row["Negative MHC alleles"] = row["Negative MHC alleles"]
        flat_row["Total response T cell assay(s)"] = row[
            "Total response T cell assay(s)"
        ]

        # Add one column per T cell assay type, like the old script did.
        for assay_name in assay_column_names:
            flat_row[assay_name] = row["assay_entries"].get(assay_name, "-")

        flat_row["Source"] = row["Source"]
        flat_rows.append(flat_row)

    return flat_rows


def make_long_rows(scraped_rows, assay_column_names, organism_override=None):
    """
    Convert scraped rows into a LONG table.

    Why long?
    - Best for analysis in pandas or R.
    - One row = one antigen-organism pair.
    - Easier to filter, group, and count later.

    In the long table, the MHC/T cell summary is repeated for each pair.
    That is normal, because those summaries belong to the epitope page.
    """
    long_rows = []

    for row in scraped_rows:
        # If there are no pairs, keep one placeholder row so we do not lose the source.
        if not row["pairs"]:
            long_row = OrderedDict(
                {
                    "Epitope": row["Epitope"],
                    "Pair_index": "-",
                    "Antigen": "-",
                    "UniProt": "-",
                    "Organism": "-",
                    "Positive MHC alleles": row["Positive MHC alleles"],
                    "Negative MHC alleles": row["Negative MHC alleles"],
                    "Total response T cell assay(s)": row[
                        "Total response T cell assay(s)"
                    ],
                }
            )

            for assay_name in assay_column_names:
                long_row[assay_name] = row["assay_entries"].get(assay_name, "-")

            long_row["Source"] = row["Source"]
            long_rows.append(long_row)
            continue

        # Otherwise create one row per pair.
        for i, pair in enumerate(row["pairs"], start=1):
            # long_row = OrderedDict(
            #     {
            #         "Epitope": row["Epitope"],
            #         "Pair_index": i,
            #         "Antigen": pair["antigen"],
            #         "UniProt": pair["uniprot"],
            #         "Organism": organism_override or pair["organism"],
            #         "Positive MHC alleles": row["Positive MHC alleles"],
            #         "Negative MHC alleles": row["Negative MHC alleles"],
            #         "Total response T cell assay(s)": row["Total response T cell assay(s)"],
            #     }
            # )
            long_row = OrderedDict(
                {
                    "Epitope": row["Epitope"],
                    "Pair_index": i,
                    # These values come from the antigen-organism pair.
                    # .get("key", "-") means:
                    # "Use this value if it exists; otherwise use '-'."
                    # i.e. Using .get() so the script does not crash if one key is missing.
                    # Example: if a pair has no UniProt code, we write "-" instead.
                    # They are repeated for each antigen-organism pair in the long table.
                    "Antigen": pair.get("antigen", "-"),
                    "UniProt": pair.get("uniprot", "-"),
                    "Organism": organism_override or pair.get("organism", "-"),
                    # These columns belong to the epitope page, so we repeat them for each pair.
                    "Positive MHC alleles": row["Positive MHC alleles"],
                    "Negative MHC alleles": row["Negative MHC alleles"],
                    "Total response T cell assay(s)": row[
                        "Total response T cell assay(s)"
                    ],
                }
            )

            # Add the dynamic T-cell assay columns.
            # Example: "qualitative binding" -> "1/1"
            for assay_name in assay_column_names:
                long_row[assay_name] = row["assay_entries"].get(assay_name, "-")

            # Keep the source URL at the end.
            long_row["Source"] = row["Source"]

            # Add this completed row to the long table.
            long_rows.append(long_row)

    return long_rows


def save_failed_link(
    failed_csv,
    source_link,
    error_message,
    failure_type="unknown",
    reference_string="-",
):
    """
    Save failed links into a small CSV file.

    Why do this?
    - The log file is useful for humans.
    - A CSV file is easier to reuse as input for a retry run.

    We store two main failure types:
    - download: the page could not be downloaded, for example 429 Too Many Requests.
    - parsing: the page was downloaded, but the script could not extract the data.
    """

    failed_csv = Path(failed_csv)

    # Make sure the output folder exists before writing the file.
    # If the folder already exists, this does nothing.
    failed_csv.parent.mkdir(parents=True, exist_ok=True)

    # Check whether the file already exists.
    # If it does not exist, we write the header first.
    file_exists = failed_csv.exists()

    with open(failed_csv, mode="a", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "Source",
                "Failure_type",
                "Error",
                "Reference_string",
            ],
        )

        # Write the header only once.
        if not file_exists:
            writer.writeheader()

        writer.writerow(
            {
                "Source": source_link,
                "Failure_type": failure_type,
                "Error": error_message,
                "Reference_string": reference_string,
            }
        )


# ---------------------------------- Main -------------------------------------


def main():
    """
    Main function:
    1. Read arguments
    2. Set up logging
    3. Read links
    4. Download and parse each page
    5. Build the wide DataFrame
    6. Build the long DataFrame
    7. Save both CSV files
    """
    args = parse_args()

    # Decide where to save the long output table.
    output_long_csv = get_long_output_path(args.output_csv, args.output_long_csv)

    # Decide where to save the failed-links CSV.
    # This puts it in the same directory as the main wide output CSV.
    failed_csv = get_failed_output_path(args.output_csv)

    # Start each run with a clean failed-links CSV.
    #
    # Why?
    # - The wide and long output CSV files are overwritten each run.
    # - So the failed-links CSV should also represent only the current run.
    failed_csv_path = Path(failed_csv)
    if failed_csv_path.exists():
        failed_csv_path.unlink()

    # --------------------------- Logging setup ---------------------------
    logging.basicConfig(
        level=logging.DEBUG,
        filename=args.log,
        filemode="w",
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%d/%m/%Y, %H:%M:%S",
    )

    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)

    # --------------------------- Read input links ------------------------
    links = get_iedb_link_list(args.input)
    logging.info(f"There are {len(links)} links to process.")

    scraped_rows = []

    # --------------------------- Process links ---------------------------
    # Build one HTTPX client and reuse it for all requests.
    with build_httpx_client() as client:
        with alive_bar(len(links)) as bar:
            for index, link in enumerate(links, start=1):
                normalized_link = normalize_iedb_url(link)

                try:
                    html = extract_html(
                        normalized_link,
                        client,
                        retry_base_delay=args.retry_base_delay,
                        max_retries=args.max_retries,
                    )
                    row = extract_all_data_from_html(
                        html,
                        normalized_link,
                        failed_csv,
                    )
                except Exception as e:
                    logging.exception(f"Failed to process link {link}: {e}")

                    # Save download-level failures too.
                    #
                    # Example download failures:
                    # - 429 Too Many Requests
                    # - timeout
                    # - connection error
                    # - HTTP 500 server error
                    save_failed_link(
                        failed_csv=failed_csv,
                        source_link=normalized_link,
                        error_message=str(e),
                        failure_type="download",
                        reference_string="-",
                    )

                    row = {
                        "Epitope": "-",
                        "pairs": [],
                        "Positive MHC alleles": "-",
                        "Negative MHC alleles": "-",
                        "Total response T cell assay(s)": "-",
                        "assay_entries": OrderedDict(),
                        "Source": normalized_link,
                    }

                scraped_rows.append(row)
                logging.info(f"Processed link {index}/{len(links)}")
                bar()

                # ------------------------------------------------------------
                # Optional normal pause between DIFFERENT links.
                # ------------------------------------------------------------
                # This is controlled by:
                #   -s / --sleep_seconds
                #
                # Example:
                #   -s 2
                # means:
                #   after finishing this link, wait 2 seconds before the next link.
                #
                # Use -s 0 to go immediately from one successful link to the next.
                #
                # We do not sleep after the last link because there is no next link.
                if args.sleep_seconds > 0 and index < len(links):
                    logging.info(
                        f"Sleeping {args.sleep_seconds} seconds before the next link."
                    )
                    time.sleep(args.sleep_seconds)

    logging.info("Finished extracting data from all links.")
    logging.info(f"Failed-links CSV saved to: {failed_csv}")

    # --------------------------- Build wide table ------------------------
    assay_column_names = collect_all_assay_column_names(scraped_rows)
    flat_rows = flatten_rows(
        scraped_rows, assay_column_names, organism_override=args.organism
    )
    wide_df = pd.DataFrame(flat_rows)
    wide_df.fillna("-", inplace=True)

    if "Antigen_1" in wide_df.columns:
        wide_df = wide_df.sort_values(by=["Antigen_1", "Epitope"]).reset_index(
            drop=True
        )
    else:
        wide_df = wide_df.sort_values(by="Epitope").reset_index(drop=True)

    # --------------------------- Build long table ------------------------
    long_rows = make_long_rows(
        scraped_rows, assay_column_names, organism_override=args.organism
    )
    long_df = pd.DataFrame(long_rows)
    long_df.fillna("-", inplace=True)
    long_df = long_df.sort_values(by=["Epitope", "Pair_index", "Antigen"]).reset_index(
        drop=True
    )

    # --------------------------- Save CSV files -------------------------
    wide_df.to_csv(args.output_csv, index=False)
    long_df.to_csv(output_long_csv, index=False)

    logging.info(f"Wide CSV saved to: {args.output_csv}")
    logging.info(f"Long CSV saved to: {output_long_csv}")


# This line means:
# "Only run main() if this file is executed directly."
if __name__ == "__main__":
    main()
