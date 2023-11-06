"""
get sequence length for each gene in the gene catalogue
"""

from Bio import SeqIO
import sys
import csv

input_f = sys.argv[1]
output_f = sys.argv[2]

gene_catalogue = SeqIO.parse(input_f, "fasta")

with open(output_f, 'w') as f:
    csv_writer = csv.writer(f)
    csv_writer.writerow(["gene", "length"])
    for record in gene_catalogue:
        csv_writer.writerow([record.id, len(record.seq)])
