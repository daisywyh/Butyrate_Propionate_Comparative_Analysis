from Bio import SeqIO
import sys

# file = "/Users/rebeccachristensen/Desktop/Cremer_Lab_2022/propionate_gene_catalogue/workflow/out/compiled_gene_catalogue.fa"
# new_file = "../out/compiled_gene_catalogue_2.fa"

file = sys.argv[1]
new_file = sys.argv[2]

records = SeqIO.parse(file, 'fasta')
new_records = []

for record in records:
    strain_gene = str(record.description).split(" ")[1].split(" ")[0]
    record.id = strain_gene
    # print(record.id)
    new_records.append(record)

# print(new_records)

SeqIO.write(new_records, new_file, "fasta")
