#!/usr/bin/env python3

# *----------------------------------------------------------------------------
# | PROGRAM NAME: iedb_scraper.py
# | DATE: 08/01/24
# | CREATED BY: Lila Maciel Rodriguez Perez
# | PROJECT FILE: IEDB Scraper
# | VERSION: 3.0
# *----------------------------------------------------------------------------

# *---------------------------------  Libraries -------------------------------

import argparse
import json
import logging
import re
from collections import OrderedDict

import pandas as pd
import requests
from alive_progress import alive_bar
from bs4 import BeautifulSoup


def parse_args():
    parser = argparse.ArgumentParser(description="Parse HTML text and generate table.")
    parser.add_argument("-i", "--input", required=True, help="File with links.")
    parser.add_argument("-o", "--output_csv", required=True, help="Output CSV file.")
    parser.add_argument("-org", "--organism", help="Specify the organism.")
    parser.add_argument(
        "-l",
        "--log",
        default="log.txt",
        help='Log file name (optional). Default is "log.txt".',
    )
    return parser.parse_args()


def get_iedb_link_list(file):
    links_list = []
    with open(file, mode="r") as fh:
        for link in fh:
            link = link.strip()
            links_list.append(link)
    return links_list


def extract_html(url):
    response = requests.get(url)
    return response.text


def get_epitope_data(pattern, string):
    # Using regular expressions to find the JSON-like data
    match = re.search(pattern, string)  # epitope_data_
    # If there is no match
    if not match:
        raise ValueError("There was no match for epitope data in the HTML")
    # Get epitope data and Convert the matched strings to valid JSONs
    epitope_data_string = json.loads(match.group(1))["data"]["referenceEpitopeString"]
    epitope_match = re.search(
        r"(^[A-Z]+(?: \+ [A-Z]+\([A-Z0-9]+\))?)\s.*", epitope_data_string
    )
    antigen_match = re.search(
        r".*studied as part of (.*?) from .*", epitope_data_string
    )
    organism_match = re.search(r".*from (.*?)(?:\.|$)", epitope_data_string)
    # If there is not a match for the tree of them:
    if not (epitope_match and antigen_match and organism_match):
        raise ValueError("There was not a match for epitope, antigen, or organism")
    # If there was a match for the tree elements, return them as a tuple
    return organism_match.group(1), antigen_match.group(1), epitope_match.group(1)


def get_json_data(pattern, string):
    # Using regular expressions to find the JSON-like data
    compiled_data_match = re.search(pattern, string)
    # If there was no match
    if not compiled_data_match:
        raise ValueError("There was no match in the HTML")
    # If there was a match, proceed
    try:
        return json.loads(compiled_data_match.group(1).replace("'", '"'))
    except json.JSONDecodeError:
        raise ValueError("Failed to parse data as JSON")


def get_allele_data(json_data):
    """
    Extracts MHC molecule data from a given JSON structure.

    Args:
    json_data (dict): A dictionary representing JSON data.

    Returns:
    list: A list of dictionaries with MHC molecule information if available, otherwise an empty list.
    """
    try:
        # Ensure that the keys 'data' and the nested 'data' exist
        if "data" in json_data and "data" in json_data["data"][0]:
            # List of dictionaries with the MHC molecules in each dictionary
            return json_data["data"][0]["data"]
        else:
            return []
    except (KeyError, IndexError, TypeError):
        # Handle the case where the structure is not as expected
        return []  # or raise an appropriate exception


def allele_clasification(mhc_list):
    """
    Classifies MHC molecules based on their positive count into positives and negatives.

    Args:
    mhc_list (list of dict): List of dictionaries with the MHC molecules in each dictionary.

    Returns:
    tuple: Two lists containing positive and negative MHC molecules.
    """
    # mhc_list: List of dictionaries with the MHC molecules in each dictionary
    # List variables
    mhc_positives = []
    mhc_negatives = []
    # Iterating over each dictionary-element of the list
    for mhc_dict in mhc_list:
        if "mhc_molecule" in mhc_dict and "positive_count" in mhc_dict:
            try:
                mhc_molecule = mhc_dict["mhc_molecule"]
                mhc_positive_count = int(mhc_dict["positive_count"])
                # Check the positive count
                if mhc_positive_count > 0:
                    mhc_positives.append(mhc_molecule)
                else:
                    mhc_negatives.append(mhc_molecule)
            except ValueError:
                # Handle the case where positive_count is not an integer
                continue  # Skip this entry if positive_count is not a valid integer
        else:
            continue  # Skip this entry if required keys are missing
    return mhc_positives, mhc_negatives


def get_T_cell_assay_data(json_data):
    # json_data["data"]?
    # List of dictionaries, I only want the last dict, whose 'data' key has the info I want:
    try:
        if "data" in json_data and "data" in json_data["data"][-1]:
            # list of dictionaries with T cell assays info
            return json_data["data"][-1]["data"]
        else:
            return []
    except (KeyError, IndexError, TypeError):
        # Handle the case where the structure is not as expected
        return []  # or raise an appropriate exception


def calculate_total_T_response(assay_list):
    """Calculates if the total T cell response is positive (1) or negative (0).

    Args:
        assay_list (list): List of tuples (assay type, assay result)

    Returns:
        integer: Total response of all the T cell assays
    """
    if not assay_list:
        return "-"
    # Generate list of booleans for each assay to determine if there's any positive response
    t_assays_boolean_list = [int(val.split("/")[0]) > 0 for key, val in assay_list]
    # Evaluate the list of booleans and return total response
    return 1 if any(t_assays_boolean_list) else 0


def main():
    args = parse_args()
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Get the log file name from parsed arguments
    log_file_name = args.log
    # Initialize logging (To log only to a file)
    logging.basicConfig(
        level=logging.DEBUG,
        filename=log_file_name,
        filemode="w",
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    # Create a console handler and set level to INFO
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # Set a format which is simpler for console use
    time_format = "%d/%m/%Y, %H:%M:%S"  # custom format
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt=time_format
    )
    # Tell the console handler to use this format
    console.setFormatter(formatter)
    # Add the console handler to the root logger
    logging.getLogger("").addHandler(console)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Reading the file and get the list of links
    links = get_iedb_link_list(file=args.input)

    print()
    logging.info(f"There are {len(links)} links to extract information from ...\n")

    # List to store OrderedDicts for each link
    ordered_dicts = []

    # (1) Iterate through each link in links, extract information and store it in dictionary
    with alive_bar(len(links)) as bar:
        for index, link in enumerate(links):
            # time.sleep(1)  # Sleep for 1 seconds
            html = extract_html(url=link)

            # To parse the html
            soup = BeautifulSoup(html, "html.parser")
            scripts = soup.find_all("script", type="text/javascript")

            # Dictionary variables
            epitope_dict = OrderedDict()
            # Flags
            data_extracted = False  # Flag to indicate successful data extraction

            for script in scripts:
                # Making sure that the script tag contains a string
                # and also contains the strings with the data I want
                if (
                    script.string
                    and "refernceEpitopeData" in script.string
                    and "compiledData" in script.string
                ):
                    try:
                        # Part 1: Get epitope data (organism, antigen, epitope)
                        organism, antigen, epitope = get_epitope_data(
                            r"var refernceEpitopeData = (.*?});", script.string
                        )
                        # # Part 2.0: Search the data and convert it to json format
                        compiled_data_json = get_json_data(
                            r"var compiledData = (.*?});", script.string
                        )
                        # Part 2: Get MHC allele data (positive and negative list)
                        mhc_list = get_allele_data(compiled_data_json)
                        positive_list, negative_list = allele_clasification(mhc_list)
                        # Part 3: Get T cell assay data
                        assay_list = get_T_cell_assay_data(compiled_data_json)
                        assay_info_list = [
                            (
                                assay_dict["assay_type"],
                                f"{assay_dict['positive_count']}/{assay_dict['total_count']}",
                            )
                            for assay_dict in assay_list
                        ]
                        total_T_response = calculate_total_T_response(assay_info_list)
                        data_extracted = True
                        break
                        # out of the loop, data was found, do not keep iterating over other scripts anymore
                    except ValueError as e:
                        # Handle the ValueError but continue the loop
                        logging.exception(f"> A ValueError occurred: {e}")
                    except Exception as e:
                        # Handle the Error but continue the loop
                        logging.exception(f"> An unexpected error occurred: {e}")

            if not data_extracted:
                logging.info(f"No data extracted from the link: {link}")

            if data_extracted:
                # Adding Epitope data to dictionary
                epitope_dict["Organism"] = organism
                epitope_dict["Antigen"] = antigen
                epitope_dict["Epitope"] = epitope
                # Adding MHC allele data to dictionary
                epitope_dict["Positive MHC alleles"] = (
                    ",".join(positive_list) if positive_list else "-"
                )
                epitope_dict["Negative MHC alleles"] = (
                    ",".join(negative_list) if negative_list else "-"
                )
                # Adding T Cell assays data to dictionary
                epitope_dict["Total response T cell assay(s)"] = total_T_response
                assay_entries = {
                    assay_key: assay_val
                    for assay_key, assay_val in assay_info_list
                    if assay_key not in epitope_dict
                }
                epitope_dict.update(assay_entries)

                # Adding the link to verify information
                epitope_dict["Source"] = link

            # Appending dictionary with information from this html to list of dictionaries
            ordered_dicts.append(epitope_dict)

            # Optionally, let's log a custom message with each update
            logging.info(f"Processed link {index + 1}/{len(links)}\n")

            # After processing each link, update the progress bar
            bar()

    print()
    logging.info("The data has been extracted from all the links.\n")

    print()
    # (2) Create a DataFrame from the list of OrderedDicts
    logging.info("Creating the pandas DataFrame.\n")
    # Creating df without the 'Source' column
    df = pd.DataFrame(ordered_dicts).drop(columns="Source")
    # Adding the 'Source' column to the df. It will appear as the last column.
    df["Source"] = [d["Source"] for d in ordered_dicts]
    # If an organism was specified
    if args.organism:
        df["Organism"] = args.organism
    # Replacing NaN with "-"
    df.fillna("-", inplace=True)
    # Sorting Dataframe by 'Antigen' and resseting the index
    sorted_by_antigen_df = df.sort_values(by="Antigen").reset_index(drop=True)
    print()

    # (3) Storing the pandas DataFrame into a CSV file
    logging.info("Storing the pandas DataFrame into a CSV file.\n")
    sorted_by_antigen_df.to_csv(args.output_csv)


if __name__ == "__main__":
    main()
