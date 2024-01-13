# IEDB Scraper

## Description

This tool is designed to automate the extraction of immunological data from the Immune Epitope Database (IEDB). It parses web pages corresponding to given IEDB links and extracts relevant information such as the organism, antigen, epitope, MHC alleles, and T cell assays. The data is then structured into a pandas DataFrame and saved as a CSV file for further analysis or reporting.

## Installation

To run this script, you will need Python installed on your system. The script uses several Python libraries, which you can install using the following commands:

```bash
pip install pandas
pip install requests
pip install bs4
pip install alive-progress
```

Ensure you have these libraries installed before executing the script.


## Usage

To use this script, you need to have a text file containing IEDB links, one per line. You can then run the script with the following command:

```bash
python script_name.py -i input_file.txt -o output_data.csv -l log_file.txt
```

- `-i` or `--input`: (Required) Path to the input file containing IEDB links.
- `-o` or `--output_csv`: (Required) Path to the output CSV file where the data will be saved.
- `-l` or `--log`: (Optional) Path to the log file. Defaults to "log.txt" if not specified.

## Example DataFrame

Below is an example of the DataFrame structure you can expect as output. This table will be saved in CSV format:

| Organism   | Antigen             | Epitope            | Positive MHC alleles         | Negative MHC alleles         | Total response T cell assay(s) | Qualitative binding | T cell binding | IFNg release | Source                               |
|------------|---------------------|--------------------|-----------------------------|-----------------------------|--------------------------------|---------------------|----------------|--------------|--------------------------------------|
| SARS-CoV1  | Replicase polyprotein 1ab | AAISDYDYY          | HLA-A*03:01, HLA-A*11:01    | HLA-A*31:01, HLA-A*33:01    | 1                              | 1/1                 | NaN            | NaN          | https://www.iedb.org/epitope/234     |
| SARS-CoV2  | Nucleoprotein       | AEGSRGGSQA         | HLA-B*45:01                 | HLA-B*18:01, HLA-B*40:01    | 1                              | NaN                 | 3/3            | NaN          | https://www.iedb.org/epitope/956     |
| ...        | ...                 | ...                | ...                         | ...                         | ...                            |



*Note: The actual output will contain more rows and columns based on the extracted data.*

### Contributing

Contributions to this project are welcome. Please ensure you follow the guidelines outlined in our CONTRIBUTING.md file.

### License

This project is released under the MIT License. See the LICENSE file for more details.

### Contact

For any queries or discussions regarding this tool, please open an issue in the repository, and we will

get back to you promptly.

### Acknowledgments

This tool was created to assist researchers and professionals in the field of immunology. A big thank you to the Immune Epitope Database and Analysis Resource (IEDB) for maintaining the database and providing a valuable resource to the scientific community.