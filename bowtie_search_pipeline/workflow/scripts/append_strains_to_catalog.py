strains_f = "config/housekeeping_strains_add_to_catalog.txt"
compiled_genomes_f = "workflow/out/compiled_nuc_prodigal_translation_hit_strains_appended.fna"

with open(strains_f, 'r') as f:
    strains = f.read()
    strains = strains.split()

for strain in strains:
    genome_f = "workflow/out/nucProdigalwHeader/" + strain + "_prodigal_nuc_addHeader.fna"

    compiled_genomes_file = open(compiled_genomes_f, 'a+')
    genome_file = open(genome_f, 'r')

    compiled_genomes_file.write(genome_file.read())

    compiled_genomes_file.close()
    genome_file.close()
