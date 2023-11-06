import pandas as pd

file = "workflow/out/pathway_abundance_calculation/bt_hits_long_withGeneLengthCorr_noDupRows.csv"
output_file = "workflow/out/pathway_abundance_calculation/pathway_abundance_calculation_batch_commands.txt"

df = pd.read_csv(file)

n_rows = len(df.index)
n_rows_per_batch = 50000
n_batches = n_rows // n_rows_per_batch


command = "\npython3 workflow/scripts/pathway_abundance_calculation_tuples.py "

commands_list = []

for i in range(0, n_rows, n_rows_per_batch):
    # check for remainder
    if n_rows - i < n_rows_per_batch:
        batch_command = command + " " + str(i) + " " + str(n_rows)

        commands_list.append(batch_command)
    else:
        index_start = i
        index_end = i + n_rows_per_batch

        batch_command = command + " " + str(index_start) + " " + str(index_end) + " &"

        commands_list.append(batch_command)

f = open(output_file, 'w')
f.writelines(commands_list)
f.close()
