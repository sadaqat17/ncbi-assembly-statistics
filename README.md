
# NCBI Assembly Statistics Fetcher

This repository contains a Python script to fetch genome assembly statistics for multiple species or genera using NCBI's Entrez API.

## Features

- Fetch statistics for all genome assemblies associated with a given TaxID or species name.
- Includes detailed statistics such as genome size, assembly level, N50 values, sequencing technology, and more.
- Outputs the results to a CSV file.

## Prerequisites

- Python 3.x
- Biopython library (for NCBI Entrez API interactions)
- An internet connection to access NCBI's databases

## Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/sadaqat17/ncbi-assembly-statistics.git
   cd ncbi-assembly-statistics
   ```

2. Install the required Python library:

   ```bash
   pip install biopython
   ```

3. Ensure that the `fetch_assemblies.py` script and the `species_list.txt` file are in the repository folder.

## Usage

### Prepare the Species List

Open `species_list.txt` and add the species names or TaxIDs you want to analyze. Each species should be on a new line. The script supports both species names and TaxIDs. Example `species_list.txt`:

```
Homo sapiens
3702
Arabidopsis thaliana
9606
```

### Run the Script

Open a terminal, navigate to the repository folder, and run:

```bash
python fetch_assemblies.py
```

The script will fetch genome assembly statistics for each species and save the results in a CSV file named `assemblies_statistics.csv`.

## Output

The script will generate a CSV file with detailed assembly statistics, including:

- Assembly Accession
- Organism Name
- TaxID
- Submitter
- Assembly Level
- GenBank FTP Link
- RefSeq FTP Link
- Contig N50
- Scaffold N50
- Genome Size
- Chromosome Count
- Scaffold Count
- Contig Count
- Sequencing Technology
- Assembly Method
- DOI/URL

### Sample Output

A sample output row in `assemblies_statistics.csv` will look like this:

```
AssemblyAccession, OrganismName, TaxID, Submitter, AssemblyLevel, GenBankFTPLink, RefSeqFTPLink, ContigN50, ScaffoldN50, GenomeSize, ChromosomeCount, ScaffoldCount, ContigCount, SequencingTechnology, AssemblyMethod, DOI/URL
```

## License

This project is licensed under the MIT License.

## Contributions

Contributions and suggestions are welcome! Please open an issue or submit a pull request.

## Contact

If you have any questions or feedback, feel free to contact me at:  
**Email:** muhammad.sadaqat@univ-rennes.fr
