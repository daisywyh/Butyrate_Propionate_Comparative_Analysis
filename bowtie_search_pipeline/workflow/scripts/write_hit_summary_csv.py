import csv
import sys
import ast

input_f = sys.argv[1]
output_f = sys.argv[2]

# get list of genes found in hit_summaries to write header
with open(input_f, 'r') as f:
    genes_list = []
    for line in f:
        curr_dict = ast.literal_eval(line.rstrip())
        for key in curr_dict.keys():
            if key not in genes_list:
                genes_list.append(key)

with open(output_f, 'w') as f:
    csv_writer = csv.DictWriter(f, genes_list)
    csv_writer.writeheader()

    with open(input_f, 'r') as f:
        for line in f:
            curr_dict = ast.literal_eval(line.rstrip())
            csv_writer.writerow(curr_dict)
