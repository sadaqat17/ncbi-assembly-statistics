from Bio import Entrez
import csv

# Set your email for NCBI
Entrez.email = "your_email@example.com"

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

# Taxonomic ID to query
taxid = "9606"  # Example: Homo sapiens

# Fetch all assembly UIDs for the given TaxID
assembly_uids = fetch_assemblies_for_taxid(taxid)
print(f"Found {len(assembly_uids)} assemblies for TaxID {taxid}")

# Prepare output file
output_file = "assemblies_statistics.csv"
fields = [
    "Assembly Accession", "Organism Name", "TaxID", "Submitter", "Assembly Level",
    "GenBank FTP Link", "RefSeq FTP Link", "Contig N50", "Scaffold N50",
    "Base Chromosome No.", "Total Chromosomes", "Ploidy Level", "Genome Size",
    "Genome Type", "Chromosome Count", "Scaffold Count", "Contig Count",
    "Ensembl", "Phytozome", "Assembly Version", "Sequencing Technology",
    "Assembly Method", "First Reported Paper", "Year", "DOI/URL"
]

# Write statistics to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)  # Write header

    for uid in assembly_uids:
        details = fetch_assembly_details(uid)
        docsum = details['DocumentSummarySet']['DocumentSummary'][0]

        row = [
            safe_extract(docsum, "AssemblyAccession"),
            safe_extract(docsum, "Organism"),
            safe_extract(docsum, "Taxid"),
            safe_extract(docsum, "Submitter"),
            safe_extract(docsum, "AssemblyStatus"),
            safe_extract(docsum, "FtpPath_GenBank"),
            safe_extract(docsum, "FtpPath_RefSeq"),
            safe_extract(docsum, "ContigN50"),
            safe_extract(docsum, "ScaffoldN50"),
            safe_extract(docsum, "Chromosome", "Not Available"),
            safe_extract(docsum, "TotalChromosomes", "Not Available"),
            safe_extract(docsum, "Ploidy", "Not Available"),
            safe_extract(docsum, "GenomeSize", "Not Available"),
            safe_extract(docsum, "GenomeType", "Not Available"),
            safe_extract(docsum, "ChromosomeCount", "Not Available"),
            safe_extract(docsum, "ScaffoldCount", "Not Available"),
            safe_extract(docsum, "ContigCount", "Not Available"),
            safe_extract(docsum, "Ensembl", "Not Available"),
            safe_extract(docsum, "Phytozome", "Not Available"),
            safe_extract(docsum, "AssemblyVersion", "Not Available"),
            safe_extract(docsum, "SequencingTechnology", "Not Available"),
            safe_extract(docsum, "AssemblyMethod", "Not Available"),
            safe_extract(docsum, "FirstReportedPaper", "Not Available"),
            safe_extract(docsum, "Year", "Not Available"),
            safe_extract(docsum, "DOI", "Not Available"),
        ]

        print(f"Processed: {row[0]} ({row[1]})")
        writer.writerow(row)

print(f"Assembly statistics saved to {output_file}")
