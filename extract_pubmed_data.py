import pandas as pd
from Bio import Entrez
import os

# Set NCBI Entrez API credentials
Entrez.email = "arfinphd@gmail.com"  # Provide your email address (required by NCBI)
Entrez.api_key = "ce7ef4f77e3404a012ce5f96ba95b466f508"  # Your NCBI API key

# Define a function to fetch the abstract of a PubMed article using its PMID
def fetch_article_by_pmid(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        handle.close()
        return abstract
    except Exception as e:
        print(f"Error fetching PMID {pmid}: {e}")
        return None

# Read the input Excel file containing PMIDs and article titles
excel_file = "C:/dataanalysis/pythonProject1/pubmed.pigdata2.xlsx"  # Path to the Excel file
df = pd.read_excel(excel_file)
print(df.columns)  # Print column names to confirm structure

# Clean the PMIDs column by converting to string and removing decimal points
df['PMIDs'] = df['PMIDs'].astype(str).str.split('.').str[0]

# Create a directory to save downloaded article abstracts (if it doesn't exist)
output_dir = "downloaded_articles"
os.makedirs(output_dir, exist_ok=True)

# Iterate over each row in the DataFrame to fetch and save article abstracts
for index, row in df.iterrows():
    pmid = row['PMIDs']  # Get the PMID
    title = row['Title']  # Get the article title

    print(f"Fetching article for PMID: {pmid}")
    abstract = fetch_article_by_pmid(pmid)  # Fetch the abstract using Entrez

    if abstract:
        # Save the abstract and title to a text file
        filename = f"{output_dir}/PMID_{pmid}.txt"
        with open(filename, "w", encoding="utf-8") as file:
            file.write(f"Title: {title}\n\n")
            file.write(abstract)

        print(f"Saved: {filename}")
    else:
        print(f"Failed to fetch article for PMID: {pmid}")

print(f"All articles saved in '{output_dir}' directory.")
