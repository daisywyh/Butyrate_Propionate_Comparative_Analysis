rule nucProdigal:
    input:
        join(config["strainGenomeDir"],"{hit_strain}.fa")
    output:
        join(config["nucProdigalDir"],"{hit_strain}_prodigal_nuc.fna")
    conda:
        config["prodigalEnv"]
    shell:
        """
        prodigal -i 
        {input} -d {output}
        """

rule addHeaderInfo:
    input:
        join(config["nucProdigalDir"],"{hit_strain}_prodigal_nuc.fna")
    output:
        join(config["nucProdigalDirwHeader"],"{hit_strain}_prodigal_nuc_addHeader.fna")
    shell:
        """
        python3 workflow/scripts/add_fasta_header.py {input} {output} {wildcards.hit_strain}
        """

rule makeGeneCatalogue:
    input:
        strain_hits_csv_f="config/{pathway}/0.5scoreFilter_nGenes_{pathway}_HMMER_hits.csv"
    output:
        "workflow/out/{pathway}/{pathway}_{gene}_catalogue.faa"
    shell:
        """
        python3 workflow/scripts/make_gene_catalogue.py {input.strain_hits_csv_f} {output} {wildcards.gene}
        """

rule compileGeneCatalogue:
    input:
        expand("workflow/out/{pathway}/{pathway}_{gene}_catalogue.faa", pathway=PATHWAY, gene=GENES)
    output:
        "workflow/out/compiled_pathway_catalogues/compiled_{pathway}_catalogue.faa"
    shell:
        """
        cat {input} > {output}
        """

