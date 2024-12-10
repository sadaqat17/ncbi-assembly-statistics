from Bio import Entrez
import csv

# Set your email for NCBI
Entrez.email = "your_email@example.com"

def get_taxid_from_species(species):
    """
    Convert species name to TaxID.
    """
    handle = Entrez.esearch(db="taxonomy", term=species, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if record["IdList"]:
        return record["IdList"][0]  # Return the first match
    else:
        return None  # No TaxID found

def fetch_assemblies_for_taxid(taxid):
    """
    Fetch all genome assembly UIDs for a given TaxID.
    """
    handle = Entrez.esearch(db="assembly", term=f"txid{taxid}[Organism:exp]", retmax=1000)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_assembly_details(uid):
    """
    Fetch detailed assembly information for a given UID.
    """
    handle = Entrez.esummary(db="assembly", id=uid, report="full")
    record = Entrez.read(handle)
    handle.close()
    return record

def safe_extract(docsum, key, default="Not Available"):
    """
    Safely extract a value from the document summary or return a default.
    """
    return docsum.get(key, default)

# Read species from the species_list.txt file
with open("species_list.txt", "r") as file:
    species_list = file.readlines()

# Clean up the species names (remove newlines and extra spaces)
species_list = [species.strip() for species in species_list]

# Prepare output file
output_file = "assemblies_statistics.csv"
fields = [
    "Assembly Accession", "Organism Name", "TaxID", "Assembly Level",
    "GenBank FTP Link", "RefSeq FTP Link", "Contig N50", "Scaffold N50"
]

# Write statistics to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)  # Write header

    # Loop through the species list
    for species in species_list:
        print(f"Processing: {species}")

        # Check if the input is a TaxID or species name
        if species.isdigit():  # If it's a TaxID (a number)
            taxid = species
        else:  # If it's a species name
            taxid = get_taxid_from_species(species)
            if taxid is None:
                print(f"TaxID not found for {species}. Skipping.")
                continue

        print(f"Fetching data for TaxID: {taxid}")
        assembly_uids = fetch_assemblies_for_taxid(taxid)

        for uid in assembly_uids:
            details = fetch_assembly_details(uid)
            docsum = details['DocumentSummarySet']['DocumentSummary'][0]

            row = [
                safe_extract(docsum, "AssemblyAccession"),
                safe_extract(docsum, "Organism"),
                taxid,
                safe_extract(docsum, "AssemblyStatus"),
                safe_extract(docsum, "FtpPath_GenBank"),
                safe_extract(docsum, "FtpPath_RefSeq"),
                safe_extract(docsum, "ContigN50"),
                safe_extract(docsum, "ScaffoldN50")
            ]

            print(f"Processed: {row[0]} ({row[1]})")
            writer.writerow(row)

print(f"Assembly statistics saved to {output_file}")
