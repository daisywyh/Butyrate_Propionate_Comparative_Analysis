""""
so i guess i write comments now. this code is to making a goddamn gene catalogue from a given fasta file of compiled
HMMER hits and a csv file of strains that passed the filtering process.

1. iterate through the list of hit strains in the csv_file
2. open each strain fasta file
3. iterate through the strain fasta file records and match for the sequence description
4. add the sequence that matched in (3) and append to new set of records for the gene catalog
5. write new records in (4)
"""

import csv
import sys
from Bio import SeqIO

# input
csv_file = sys.argv[1]
# output
output_fasta = sys.argv[2]
# param
gene = sys.argv[3]

gene_catalogue_records = []
sequence_catalogue = []

with open(csv_file, 'r') as f:
    dict_reader = csv.DictReader(f)
    for row in dict_reader:
        row_description = row[("hit_id." + gene)] + " " + row[("description." + gene)]
        row_gene_coord = row_description.split("/")[0]
        row_strain = row["strain"]
        #print("row strain: " + row_strain)

        genome_file = "workflow/out/nucProdigalwHeader/" + row_strain + "_prodigal_nuc_addHeader.fna"
        for record in SeqIO.parse(genome_file, "fasta"):
            #print("genome file opened for: " + row_strain)
            record_gene_coord = (str(record.id)).split(";")[0]
            if row_gene_coord == record_gene_coord:
                #print("gene coord matched for: " + row_strain)
                coord = row_description.split("/")[1].split(" ")[0]
                start = int(coord.split("-")[0]) - 1
                end = int(coord.split("-")[1]) - 1
                record.id = row_description
                record.seq = record.seq[start:end]

                if str(record.seq) not in sequence_catalogue:
                    #print("sequence added to catalogue for: " + row_strain)
                    sequence_catalogue.append(str(record.seq))
                    gene_catalogue_records.append(record)

SeqIO.write(gene_catalogue_records, output_fasta, "fasta")
