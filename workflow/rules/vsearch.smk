###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################


rule vsearch_merge_readpairs:
    # merge read pairs after QC + trimming
    input:
        read1="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq.gz",
        read2="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[1] + ".fastq.gz",
    output:
        merged="results/02_MERGED/{sample}_merged.fastq",
    resources:
        mem_mb=8000,
    conda:
        "../envs/vsearch.yaml"
    params:
        relabel="{sample}.",
    threads: 4
    log:
        "logs/vsearch_merge/{sample}.log",
    shell:
        """
        vsearch --fastq_mergepairs {input.read1} --reverse {input.read2} --fastqout {output.merged} --relabel {params.relabel} 2> {log}
        """

rule concatenate_merged_readpairs:
    # pool samples together in one .fastq
    input:
        expand("results/02_MERGED/{sample}_merged.fastq", sample=SAMPLES.index),
    output:
        "results/02_MERGED/all_merged.fastq",
    log:
        "logs/vsearch_pool_all/pooling.log",
    shell:
        """
        cat {input} >> {output} 2> {log}
        """

rule vsearch_quality_filtering:
    # strip primers based on known length, QC based on maxee
    input:
         "results/02_MERGED/all_merged.fastq",
    output:
         "results/03_STRIPPED_FILTERED/stripped_filtered.fasta",
    resources:
        mem_mb=8000,
    conda:
        "../envs/vsearch.yaml"
    threads: 4
    log:
        "logs/vsearch_stripped_filtered/stripping_filtering.log",
    shell:
        """
        vsearch --fastq_filter {input} --fastq_stripleft 20 --fastq_stripright 17 --fastq_maxee 1 --fastq_maxlen 400 --fastq_minlen 300 --fastaout {output} 2> {log} 
        """

rule vsearch_unique:
    # dereplication, identify unique sequences
    input:
        "results/03_STRIPPED_FILTERED/stripped_filtered.fasta",
    output:
        "results/04_UNIQUE/unique.fasta",
    resources:
         mem_mb=8000,
    conda:
        "../envs/vsearch.yaml"
    threads: 4
    log:
        "logs/vsearch_unique/unique.log",
    shell:
        """
        vsearch --derep_fulllength {input} -sizeout -relabel Uniq -output {output} 2> {log}
        """

rule vsearch_otu:
    # otu clustering and chimera removal
    input:
        "results/04_UNIQUE/unique.fasta",
    output:
        "results/05_OTU_CLUSTERING/otus.fasta",
    resources:
        mem_mb=8000,
    conda:
        "../envs/vsearch.yaml"
    threads: 4
    log:
        "logs/vsearch_otus/otus.log",
    shell:
        """
        vsearch --cluster_size {input} --threads {threads} --id 0.8 -relabel OTU_ --centroids {output} 2> {log}
        """

rule vsearch_otu_table:
    # generate otu table
    input:
        stripped="results/03_STRIPPED_FILTERED/stripped_filtered.fasta",
        otus="results/05_OTU_CLUSTERING/otus.fasta",
    output:
        "results/05_OTU_CLUSTERING/otutab.txt",
    resources:
        mem_mb=8000,
    conda:
        "../envs/vsearch.yaml"
    threads: 4
    log:
        "logs/vsearch_otus/otutab.log",
    shell:
        """
        vsearch -usearch_global {input.stripped} --db {input.otus} --id 0.9 --otutabout {output} 2> {log}
        """

rule vsearch_tax_summary:
    # assign taxonomoy 
    input:
        otus="results/05_OTU_CLUSTERING/otus.fasta",
        otutab="results/05_OTU_CLUSTERING/otutab.txt",
    output:
        sintax="results/06_TAXONOMY/sintax.txt",
    resources:
        mem_mb=8000,
    conda:
        "../envs/vsearch.yaml"
    params:
        ref_db=REF_DB,
    threads: 4
    log:
        "logs/vsearch_taxonomy/sintax.log",
    shell:
        """
        vsearch -sintax {input.otus} -db {params.ref_db} -strand both -sintax_cutoff 0.5 -tabbedout {output.sintax} 2> {log}
        sed -i.tmp 's/#OTU ID//' {input.otutab}
        """

rule clean_up_taxonomy:
    # cleanup taxonomy 
    input:
        raw_tax="results/06_TAXONOMY/sintax.txt",
    output:
        clean_tax="results/06_TAXONOMY/sintax_clean.txt",
    log:
        "logs/vsearch_taxonomy/clean_up.log",
    shell:
        """
        workflow/scripts/convert_usearch_tax.sh {input.raw_tax} {output.clean_tax} 2> {log}
        """