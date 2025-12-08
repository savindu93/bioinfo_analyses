"""
Input - Genbank accessions of microarray probe IDs

This python script carries out the following,
- Use genbank accessions of probe IDs of microarray chips to search for their corresponding gene (gene ID) using NCBIWWW of Biopython package
- The gene IDs are then mapped to their corresponding entrez IDs with entrez esearch

Output - Dataframe of gene IDs with their Entrez IDs
"""



import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
import re

deg_na = pd.read_csv('deg_na.csv', index_col = 0)
# print(deg_na)


gene_info = {
    "acc": [],
    "gene_id": []
}


for acc in deg_na["ID"]:

    # Retrieve XML blastn entries of genes using their GB accessions (given as ID) of probe IDs in the deg_na.csv file
    res_handle = NCBIWWW.qblast(
        program = "blastn",
        database = "refseq_rna",
        sequence = acc,
        hitlist_size = 10,
        expect = 0.001,
        format_type = "XML",
        entrez_query = "human[organism]"
    )

    # Save blastn results in .xml files
    with open(f"{acc}_res.xml", "w") as out:
        out.write(res_handle.read())

    res_handle.close()

    # Using the above blastn results identify the gene IDs corresponding to that accession
    # from the best hit which is the 1st alignment
    with open(f"{acc}_res.xml") as res_handle:
        blast_rec = NCBIXML.read(res_handle)
        if len(blast_rec.alignments) == 0:
            gene_info["acc"].append(acc)
            gene_info["gene_id"].append("NA")
            continue

        entry = blast_rec.alignments[0].title
        print(">", entry)

        if len(entry.split("|")) > 1:
            gene_id = entry.split("|")[4]

            # This regex pattern extracts the gene ID from the entry header
            # Ex:
            pattern = re.compile(r'\((.*?)\)')
            match = pattern.search(gene_id)
            if match:
                gene_id = match.group(1)
                print(gene_id)
            
            gene_info["acc"].append(acc)
            gene_info["gene_id"].append(gene_id)
        else:
            continue

df = pd.DataFrame(gene_info)

df.to_csv("deg_na_gene_ids.csv")


# Now map gene IDs to Entrez IDs
from Bio import Entrez
df = pd.read_csv("deg_na_gene_ids.csv")

print(df)
df["entrez_id"] = pd.NA

for index, row in df.iterrows():
    
    gene_id = row["gene_id"]
    print(gene_id)
    if pd.isna(gene_id) :
        df.at[index, "entrez_id"] = "NA"
        continue
    else:
        handle = Entrez.esearch(
            db="gene", 
            term = f"{gene_id}[sym] AND Homo sapiens[orgn]")
        record = Entrez.read(handle)
        handle.close()
        print(record)

        if record["IdList"]:
            df.at[index, "entrez_id"] = record["IdList"][0]
        else:
            df.at[index, "entrez_id"] = "NA"


df.to_csv("deg_na_gene_ids_with_entrez.csv", index=False)

