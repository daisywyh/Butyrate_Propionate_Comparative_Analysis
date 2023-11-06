rule makeMSA:
    input:
        "config/referenceGeneFiles/{pathway}/{gene}.faa"
    output:
        "config/referenceGeneMSA/{pathway}/{gene}_msa.faa"
    conda:
        config["mafftEnv"]
    shell:
        """
        mafft {input} > {output}
        """

rule buildProfileHMM:
    input:
        "config/referenceGeneMSA/{pathway}/{gene}_msa.faa"
    output:
        "config/profileHMMs/{pathway}/{gene}.HMM"
    conda:
        config["hmmerEnv"]
    shell:
        """
        hmmbuild {output} {input}
        """

# DOWNLOAD ALL STRAIN GENOMES

rule retrieveStrainGenomes:
    output:
        join(config["strainGenomeDir"], "{strain}.fa")
    conda:
        "edirect"
    shell:
        """
        esearch -db assembly -query "{wildcards.strain} [ORGN] AND representative [PROP]" \
        | elink -target nuccore -name assembly_nuccore_refseq \
        | efetch -format fasta > {output}
        """

# RUN HMMER SEARCH
#Translate to NT to AA with PRODIGAL
rule runProdigal:
    input:
        join(config["strainGenomeDir"], "{strain}.fa")
    output:
        proteinSeqs=join(config["proteinSeqDir"], "{strain}_prodigal.faa")
    conda:
        config["prodigalEnv"]
    shell:
        """
        prodigal -i {input} -a {output.proteinSeqs}
        """

# run HMMER search on all genes
rule runHMMER:
    input:
        prodigal=join(config["proteinSeqDir"], "{strain}_prodigal.faa"),
        hmm_profile="config/profileHMMs/{pathway}/{gene}.HMM"
    output:
        hmmOut=(join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.hmm.out")),
        domOut=(join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.domtblout")),
        msa=(join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.sto"))
    conda:
        config["hmmerEnv"]
    shell:
        """
        hmmsearch -o {output.hmmOut} --domtblout {output.domOut} -A {output.msa} {input.hmm_profile} {input.prodigal} 
        """

rule parseHMMER:
    input:
        (join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.domtblout"))
    output:
        (join(config["hmmOutDir"],"{pathway}/{gene}/csv_summary/{strain}_{gene}_hits.csv"))
    shell:
        """
        python3 workflow/scripts/parse_hmmer_domtable.py {input} {output}
        """

rule combineCSV:
    output:
        "workflow/out/summary/{pathway}/csv_summary/compiled_{gene}_hits.csv"
    params:
        input_dir=(join(config["hmmOutDir"],"{pathway}/{gene}/csv_summary"))
    shell:
        """
        cat {params.input_dir}/*_hits.csv > {output}
        """

# convert the HMMER output MSA from stockholm format (.sto) to fasta (.faa)
rule convertMSA_toFASTA:
    input:
        (join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.sto"))
    output:
        (join(config["hmmOutDir"],"{pathway}/{gene}/faa_summary/{strain}_{gene}_hits.faa"))
    shell:
        """
        python3 workflow/scripts/convertFASTA.py {input} {output}
        """

rule combineMSA:
    output:
        "workflow/out/summary/{pathway}/faa_summary/compiled_{gene}_hits.faa"
    params:
        input_dir=(join(config["hmmOutDir"],"{pathway}/{gene}/faa_summary"))
    shell:
        """
        cat {params.input_dir}/*.faa > {output}
        """

# COMPILE HITS
rule combineAllCSV:
    output:
        "workflow/out/summary/{pathway}/compiled_hits_{pathway}.csv"
    params:
        input_dir="workflow/out/summary/{pathway}/csv_summary"
    shell:
        """
        cat {params.input_dir}/*.csv > {output}
        """

rule compileHitInfo:
    input:
        "workflow/out/summary/{pathway}/compiled_hits_{pathway}.csv"
    output:
        "workflow/out/summary/{pathway}/maxHitScoreDF_{pathway}.csv"
    shell:
        """
        python3 workflow/scripts/compileHitInfo.py {input} {output}
        """
