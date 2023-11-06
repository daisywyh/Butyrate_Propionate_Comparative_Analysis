import sys
import pandas as pd
import json

input_file = sys.argv[1]
output_file = sys.argv[2]
pathway = sys.argv[3]
read_accession = sys.argv[4]

header = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20']
df = pd.read_csv(input_file, delim_whitespace=True, names=header)

hits_dict = {"pathway": pathway,
             "read_accession": read_accession}

gene_counts = df["3"].value_counts().to_dict()

hits_dict.update(gene_counts)

with open(output_file, 'w') as f:
    json.dump(hits_dict, f)
    f.write("\n")
