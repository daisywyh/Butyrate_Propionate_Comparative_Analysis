import pandas as pd
import sys

INDEX_START = int(sys.argv[1])
INDEX_END = int(sys.argv[2])

file = "workflow/out/pathway_abundance_calculation/bt_hits_long_withGeneLengthCorr.csv"
dups_removed_file = "workflow/out/pathway_abundance_calculation/bt_hits_long_withGeneLengthCorr_noDupRows.csv"
output_file = "workflow/out/pathway_abundance_calculation/bt_unambiguous_gene_ct_batches/bt_unambiguous_gene_cts_" + str(INDEX_START) + "-" + str(INDEX_END) + ".csv"

# do this if u wanna remove duplicate rows from df first!
df = pd.read_csv(file)
df = df.drop_duplicates()
df.to_csv(dups_removed_file, index=False)

df = pd.read_csv(dups_removed_file)

new_df = pd.DataFrame(columns=['read_accession', 'pathway', 'hits_geneLengthCorrected_unambiguous', 'gene'])

# without index batching
# for row in df.itertuples():
#     if " or " not in row[5]:
#         new_row = [row[2], row[5], row[8], row[3]]
#
#         if new_df[(new_df['read_accession'] == row[2]) &
#                   (new_df['pathway'] == row[5]) &
#                   (new_df['hits_geneLengthCorrected_unambiguous'] == row[8]) & (new_df['gene'] == row[3])].empty:
#             new_df.loc[len(new_df.index)] = new_row
#
#     elif " or " in row[5]:
#         pathways = row[5].split(" or ")
#
#         for pathway in pathways:
#             read_accession = row[2]
#
#             # subset df by read_accession and pathway
#             single_pathway_subset = df[(df["read_accession"] == read_accession) & (df["strain_pathway"] == pathway)]
#             single_pathway_sum = single_pathway_subset["hits_geneLengthCorrected"].sum()
#
#             if single_pathway_sum == 0:
#                 single_pathway_adjusted_gene_count = 0
#             else:
#                 single_pathway_adjusted_gene_count = (row[8] / single_pathway_sum)
#
#             new_row = [row[2], pathway, single_pathway_adjusted_gene_count,
#                        row[3]]
#
#             if new_df[(new_df['read_accession'] == row[2]) &
#                       (new_df['pathway'] == pathway) &
#                       (new_df['hits_geneLengthCorrected_unambiguous'] == single_pathway_adjusted_gene_count) &
#                       (new_df['gene'] == row[3])].empty:
#                 new_df.loc[len(new_df.index)] = new_row

# with index batching
#
df = df.iloc[INDEX_START:INDEX_END]
df_copy = pd.read_csv(dups_removed_file)

for row in df.itertuples():
    if " or " not in row[5]:
        new_row = [row[2], row[5], row[8], row[3]]

        if new_df[(new_df['read_accession'] == row[2]) &
                  (new_df['pathway'] == row[5]) &
                  (new_df['hits_geneLengthCorrected_unambiguous'] == row[8]) & (new_df['gene'] == row[3])].empty:
            new_df.loc[len(new_df.index)] = new_row

    elif " or " in row[5]:
        pathways = row[5].split(" or ")

        for pathway in pathways:
            read_accession = row[2]
            multiple_pathways_subset = df_copy[(df_copy["read_accession"] == read_accession) &
                                               (df_copy["strain_pathway"].isin(pathways))]
            multiple_pathways_sum = multiple_pathways_subset["hits_geneLengthCorrected"].sum()

            # subset df by read_accession and pathway
            single_pathway_subset = df_copy[(df_copy["read_accession"] == read_accession) &
                                            (df_copy["strain_pathway"] == pathway)]
            single_pathway_sum = single_pathway_subset["hits_geneLengthCorrected"].sum()

            if single_pathway_sum == 0 or multiple_pathways_sum == 0:
                single_pathway_adjusted_gene_count = 0
            else:
                single_pathway_adjusted_gene_count = (row[8] * (single_pathway_sum / multiple_pathways_sum))

            new_row = [row[2], pathway, single_pathway_adjusted_gene_count,
                       row[3]]

            if new_df[(new_df['read_accession'] == row[2]) &
                      (new_df['pathway'] == pathway) &
                      (new_df['hits_geneLengthCorrected_unambiguous'] == single_pathway_adjusted_gene_count) &
                      (new_df['gene'] == row[3])].empty:
                new_df.loc[len(new_df.index)] = new_row

                # case check with Flavonifractor_plautii_bcd
                # if row[3] == "Flavonifractor_plautii_bcd":
                #     print(row[3], pathway)
                #     print("observed gene count: " + str(row[8]))
                #     print("single pathway gene sum: " + str(single_pathway_sum))
                #     print("adjusted gene count: " + str(single_pathway_adjusted_gene_count))
                #     print("multi pathway gene sum: " + str(multiple_pathways_sum))
                #     print("denominator: " + str(single_pathway_sum / multiple_pathways_sum))

new_df.to_csv(output_file)
