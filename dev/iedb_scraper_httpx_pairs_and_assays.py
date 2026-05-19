#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# PROGRAM NAME: iedb_scraper_httpx_pairs_and_assays.py
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

    return parser.parse_args()


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


def extract_html(url, client):
    """
    Download the HTML of one webpage and return it as text.
    """
    url = normalize_iedb_url(url)
    response = client.get(url)
    response.raise_for_status()
    return response.text


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


# def extract_antigen_organism_pairs(reference_string):
#     """
#     Extract ALL antigen-organism pairs from the reference sentence.

#     We isolate the part after:
#         'studied as part of'
#     and before the next sentence:
#         '. This epitope ...'

#     Then we find repeated pieces like:
#         Antigen name (UniProt:XXXXXX) from Organism name

#     Returns a list of dictionaries like:
#         [
#             {
#                 "antigen": "Replicase polyprotein 1ab",
#                 "uniprot": "P0C6X7",
#                 "organism": "SARS-CoV1",
#                 "antigen_full": "Replicase polyprotein 1ab (UniProt:P0C6X7)"
#             }
#         ]
#     """
#     # Step 1: isolate only the part that contains the repeated pairs.
#     section_match = re.search(
#         r"studied as part of\s+(.*?)(?=\.\s+This epitope\b|$)",
#         reference_string,
#         flags=re.DOTALL,
#     )

#     if not section_match:
#         raise ValueError("Could not isolate the antigen-organism section.")

#     pairs_section = clean_whitespace(section_match.group(1))

#     # Step 2: find all repeated antigen + UniProt + organism pieces.
#     pair_pattern = re.compile(
#         r"(?:^|\s+and\s+)"  # pair starts at beginning or after the word 'and'
#         r"(?P<antigen_name>.+?)\s+"  # antigen name
#         r"\(UniProt:(?P<uniprot>[A-Za-z0-9_-]+)\)\s+from\s+"  # UniProt code
#         r"(?P<organism>.*?)"  # organism name
#         r"(?=(?:\s+and\s+.+?\s+\(UniProt:[A-Za-z0-9_-]+\)\s+from\s+)|$)",
#         flags=re.DOTALL,
#     )

#     pairs = []

#     for match in pair_pattern.finditer(pairs_section):
#         antigen_name = clean_whitespace(match.group("antigen_name"))
#         uniprot_code = clean_whitespace(match.group("uniprot"))
#         organism_name = clean_whitespace(match.group("organism"))
#         antigen_full = f"{antigen_name} (UniProt:{uniprot_code})"

#         pairs.append(
#             {
#                 "antigen": antigen_name,
#                 "uniprot": uniprot_code,
#                 "organism": organism_name,
#                 "antigen_full": antigen_full,
#             }
#         )

#     if not pairs:
#         raise ValueError("Could not extract any antigen-organism pair.")

#     return pairs


def extract_antigen_organism_pairs(reference_string):
    """
    Extract ALL antigen-organism pairs from the reference sentence.

    We isolate the part after:
        'studied as part of'
    and before the next sentence:
        '. This epitope ...'

    This version handles two cases:

    Case 1:
        Antigen name (UniProt:CODE) from Organism name

    Case 2:
        Antigen name from Organism name

    If the UniProt code is missing, we store UniProt as "-".

    Returns a list of dictionaries like:
        [
            {
                "antigen": "Replicase polyprotein 1ab",
                "uniprot": "P0C6X7",
                "organism": "SARS-CoV1",
                "antigen_full": "Replicase polyprotein 1ab (UniProt:P0C6X7)"
            }
        ]
    """
    # ------------------------------------------------------------
    # Step 1: isolate the useful part of the sentence
    # only the part that contains the repeated pairs.
    # ------------------------------------------------------------
    # Example:
    #   Full sentence:
    #       "AAISDYDYY is ... studied as part of Protein X from Virus Y.
    #        This epitope has been studied ..."
    #
    #   We only want:
    #       "Protein X from Virus Y"
    # ------------------------------------------------------------
    section_match = re.search(
        r"studied as part of\s+(.*?)(?=\.\s+This epitope\b|$)",
        reference_string,
        flags=re.DOTALL,
    )
    # If the sentence does not contain the expected section, stop here.
    if not section_match:
        raise ValueError("Could not isolate the antigen-organism section.")

    # Clean spaces/newlines so the regex is easier to apply.
    pairs_section = clean_whitespace(section_match.group(1))

    # ------------------------------------------------------------
    # Step 2: find antigen-organism pairs
    # ------------------------------------------------------------
    # find all repeated antigen + UniProt + organism pieces.
    # This regex is more flexible than the previous one.
    #
    # It can match:
    #   Replicase polyprotein 1ab (UniProt:P0C6X7) from SARS-CoV1
    #
    # It can also match:
    #   Other SARS-CoV1 protein from SARS-CoV1
    #
    # Important:
    #   The UniProt part is optional now.
    # ------------------------------------------------------------
    pair_pattern = re.compile(
        r"(?:^|\s+and\s+)"  # pair starts at beginning OR after " and "
        r"(?P<antigen_name>.+?)"  # capture antigen name
        r"(?:\s+\((?P<uniprot>UniProt:[A-Za-z0-9_-]+)\))?"  # optionally capture UniProt
        r"\s+from\s+"  # match the word "from"
        r"(?P<organism>.*?)"  # capture organism name
        r"(?=(?:\s+and\s+.+?(?:\s+\(UniProt:[A-Za-z0-9_-]+\))?\s+from\s+)|$)",  # stop before next pair or end
        flags=re.DOTALL,
    )

    pairs = []

    # ------------------------------------------------------------
    # Step 3: convert each regex match into a dictionary
    # ------------------------------------------------------------
    for match in pair_pattern.finditer(pairs_section):
        antigen_name = clean_whitespace(match.group("antigen_name"))
        organism_name = clean_whitespace(match.group("organism"))

        # The UniProt group may be missing.
        raw_uniprot = match.group("uniprot")

        if raw_uniprot:
            # Convert "UniProt:P0C6X7" into "P0C6X7"
            uniprot_code = raw_uniprot.replace("UniProt:", "")

            # Keep a combined version too, useful for human reading.
            antigen_full = f"{antigen_name} (UniProt:{uniprot_code})"
        else:
            # No UniProt code was present in the sentence.
            uniprot_code = "-"
            antigen_full = antigen_name

        pairs.append(
            {
                "antigen": antigen_name,
                "uniprot": uniprot_code,
                "organism": organism_name,
                "antigen_full": antigen_full,
            }
        )

    # ------------------------------------------------------------
    # Step 4: if nothing was found, log the exact small section
    # ------------------------------------------------------------
    # This helps debug future weird cases.
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

    Rule used in the old script:
    - if any assay has a numerator > 0, return 1
    - otherwise return 0
    - if there are no assays, return '-'
    """
    if not assay_list:
        return "-"

    t_assays_boolean_list = [int(val.split("/")[0]) > 0 for _, val in assay_list]
    return 1 if any(t_assays_boolean_list) else 0


def extract_all_data_from_html(html, source_link):
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

            # This will give me a failed_links.csv file
            save_failed_link(
                failed_csv="failed_links.csv",
                source_link=source_link,
                error_message=str(e),
                reference_string=(
                    reference_string if "reference_string" in locals() else "-"
                ),
            )
        except Exception as e:
            logging.exception(f"Unexpected error while parsing {source_link}: {e}")

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


def save_failed_link(failed_csv, source_link, error_message, reference_string="-"):
    """
    Save failed links into a small CSV file.

    This is easier to analyze than digging through log.txt.
    """

    # Check whether the file already exists.
    file_exists = Path(failed_csv).exists()

    with open(failed_csv, mode="a", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["Source", "Error", "Reference_string"],
        )

        # Write the header only once.
        if not file_exists:
            writer.writeheader()

        writer.writerow(
            {
                "Source": source_link,
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
    output_long_csv = get_long_output_path(args.output_csv, args.output_long_csv)

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
                    html = extract_html(normalized_link, client)
                    row = extract_all_data_from_html(html, normalized_link)
                except Exception as e:
                    logging.exception(f"Failed to process link {link}: {e}")
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

    logging.info("Finished extracting data from all links.")

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
