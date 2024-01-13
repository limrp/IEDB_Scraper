#!/usr/bin/env python3

# *----------------------------------------------------------------------------
# | PROGRAM NAME: iedb_scraper.py
# | DATE: 08/01/24
# | CREATED BY: Lila Maciel Rodriguez Perez
# | PROJECT FILE: IEDB Scraper
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
    parser.add_argument("-o", "--output_csv", required=True, help="Output CSV file")
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
    epitope_data_match = re.search(pattern, string)
    # If there is no match
    if not epitope_data_match:
        raise ValueError("There was no match for epitope data in the HTML")
    # Get epitope data and Convert the matched strings to valid JSONs
    epitope_data_string = json.loads(epitope_data_match.group(1))["data"][
        "referenceEpitopeString"
    ]
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


def get_and_classify_allele_data(pattern, string):
    # List variables
    mhc_positives = []
    mhc_negatives = []

    # (1) Using regular expressions to find the JSON-like MHC allele data
    compiled_data_match = re.search(pattern, string)

    # If there was no match
    if not compiled_data_match:
        raise ValueError(
            "There was no match for allele and T cell assays data in the HTML"
        )
    # If there was a match, proceed
    compiled_data_json = json.loads((compiled_data_match.group(1)).replace("'", '"'))
    # List of dictionaries with the MHC molecules in each dictionary
    mhc_list = compiled_data_json["data"][0]["data"]
    # (2) Classify the MHC molecules
    # Iterating over each dictionary-element of the list
    for mhc_dict in mhc_list:
        mhc_molecule = mhc_dict["mhc_molecule"]
        mhc_positive_count = int(mhc_dict["positive_count"])
        # Classify MHC molecules
        if mhc_positive_count > 0:
            mhc_positives.append(mhc_molecule)
        else:
            mhc_negatives.append(mhc_molecule)
    return mhc_positives, mhc_negatives  # a tuple with 2 list


def get_T_cell_assay_data(pattern, string):
    # List variables
    t_assays_list = []

    # Using regular expressions to find the JSON-like data
    compiled_data_match = re.search(pattern, string)

    # If there was no match
    if not compiled_data_match:
        raise ValueError("There was no match for T cell assays data in the HTML")
    # If there was a match, proceed
    compiled_data_json = json.loads((compiled_data_match.group(1)).replace("'", '"'))
    ## List of dictionaries, I only want the last dict, whose 'data' key has the info I want:
    ## another list of dictionaries with the T cell assay data
    assays_list = compiled_data_json["data"][-1]["data"]
    for assay_dict in assays_list:
        assay_key = assay_dict["assay_type"]
        assay_val = f"{assay_dict['positive_count']}/{assay_dict['total_count']}"
        assay_tuple = (assay_key, assay_val)
        t_assays_list.append(assay_tuple)

    # Total response of T cell assays
    ## Generate list of booleans for each assay
    t_assays_boolean_list = [int(val.split("/")[0]) > 0 for key, val in t_assays_list]
    ## Evaluate list of booleans
    if any(t_assays_boolean_list):
        # If even one element of the list is True
        total_T_response = 1
    else:
        # If none is True, all are False
        total_T_response = 0
    return t_assays_list, total_T_response


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
                        # Part 2: Get MHC allele data (positive and negative list)
                        positive_list, negative_list = get_and_classify_allele_data(
                            r"var compiledData = (.*?});", script.string
                        )
                        # Part 3: Get T cell assay data
                        t_assays_list, total_T_response = get_T_cell_assay_data(
                            r"var compiledData = (.*?});", script.string
                        )
                        # print(t_assays_list)
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
                    ", ".join(positive_list) if positive_list else "No data found"
                )
                epitope_dict["Negative MHC alleles"] = (
                    ", ".join(negative_list) if negative_list else "No data found"
                )
                # Adding T Cell assays data to dictionary
                epitope_dict["Total response T cell assay(s)"] = total_T_response
                assay_entries = {
                    assay_key: assay_val
                    for assay_key, assay_val in t_assays_list
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
    # Sorting Dataframe by 'Antigen' and resseting the index
    sorted_by_antigen_df = df.sort_values(by="Antigen").reset_index(drop=True)

    print()
    # (3) Storing the pandas DataFrame into a CSV file
    logging.info("Storing the pandas DataFrame into a CSV file.\n")
    sorted_by_antigen_df.to_csv(args.output_csv)


if __name__ == "__main__":
    main()
