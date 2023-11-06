"""
add descriptive header to nucleotide prodigal gene predictor output to include strain name
and match csv description
"""

from Bio import SeqIO
import sys

fasta_f = sys.argv[1]
output_f = sys.argv[2]
strain = sys.argv[3]

records = SeqIO.parse(fasta_f, "fasta")

records_add_header = []

for record in records:
    record.id = str(record.id) + ";" + strain
    records_add_header.append(record)

SeqIO.write(records_add_header, output_f, "fasta")

