

rule fastqc_fastq:
    input:
        PATH_FASTQ+"/{sample}.fastq.gz"
    output:
        html=PATH_QC+"/{sample}_fastqc.html",
        zip=PATH_QC+"/{sample}_fastqc.zip"
    params:
        config["fastqc"]["extra"]
    wrapper:
        "0.17.0/bio/fastqc"


rule multiqc:
    input:
        expand(PATH_QC+"/{sample}_trimmed{pe}_fastqc.zip", sample=IDs, pe=['_1','_2'])
    output:
        PATH_QC+'/multiqc.html'
    log:
        PATH_LOG+'/multiqc.log'
    wrapper:
        '0.22.0/bio/multiqc'

rule cutadapt:
    input:
        [PATH_FASTQ+"/{sample}_1.fastq.gz", PATH_FASTQ+"/{sample}_2.fastq.gz"]
    output:
        fastq1=temp(PATH_FASTQ+"/{sample}_trimmed_1.fastq.gz"),
        fastq2=temp(PATH_FASTQ+"/{sample}_trimmed_2.fastq.gz"),
        qc=PATH_QC+"/{sample}_cutadapt.txt"
    params:
        config["cutadapt"]["extra"]
    log:
        PATH_LOG + "/{sample}_cutadapt.log"
    wrapper:
        "0.23.1/bio/cutadapt/pe"

rule bwa:
    input:
        reads = [PATH_FASTQ+"/{sample}_trimmed_1.fastq.gz",PATH_FASTQ+"/{sample}_trimmed_2.fastq.gz"]
    params:
        index = config['bwa']['index'],
        sort="samtools",
        sort_order="coordinate"
    threads:
        config['Ncores']
    log:
        PATH_LOG + "/{sample}.bwa_log"
    output:
        PATH_BAM+"/{sample}.bam"
    wrapper:
        "0.23.1/bio/bwa/mem"


rule bam_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule bam_reheader:
    input:
        "{sample}.bam"
    output:
        "{sample}.rehead.bam"
    params:
        '{sample}'
    shell:
        r"""
            (samtools view -H {input}; echo -e '@RG\tID:{params}\tSM:{params}') | samtools reheader - {input}  > {output}
        """

rule bam_readgroup:
    input:
        path.join(PATH_BAM, "{sample}.bam")
    output:
        path.join(PATH_BAM, "{sample}.rg.bam")
    log:
        path.join(PATH_LOG, '{sample}_picard_rg.log')
    params:
        "RGID={sample} RGPL=illumina RGSM={sample} RGLB=library RGPU=platform"
    wrapper:
        "0.24.0/bio/picard/addorreplacereadgroups"

#rule compress_vcf:
#    input:
#        "{sample}.vcf"
#    output:
#        "{sample}.vcf.gz"
#    shell:
#        """
#            bgzip -c {input} > {output}
#        """

rule index_vcf:
    input:
        "{sample}.vcf.gz"
    output:
        "{sample}.vcf.gz.csi"
    shell:
        """
           bcftools index {input}
        """


rule merge_vcf:
    input:
        vcf=expand(PATH_VAR+"/{sample}.vcf.gz", sample=IDs),
        index=expand(PATH_VAR+"/{sample}.vcf.gz.csi", sample=IDs)

    output:
        PATH_VAR + "/merged.vcf.gz"
    shell:
        """
            bcftools merge -o {output} -O z {input.vcf:q}
        """


#rule vardict:
#    input:
#        bam=PATH_BAM + "/{sample}.bam",
#        bai=PATH_BAM + "/{sample}.bam.bai"
#    output:
#        PATH_VAR + "/{sample}.vcf"
#    log:
#        PATH_LOG + "/{sample}_vardict"
#    params:
#        bed = config['vardict']['bedfile'],
#        ref = config['vardict']['reference'],
#        path = config['vardict']['path_vardict'],
#        AF_thr = config['vardict']['allele_frequency'],
#        other_option = config['vardict']['options']
#    shell:
#        """
#          {params.path}/build/install/VarDict/bin/VarDict -G {params.ref} -f {params.AF_thr} -N {wildcards.sample} -b {input.bam} {params.other_option} {params.bed} |
#          {params.path}/VarDict/teststrandbias.R |
#          {params.path}/VarDict/var2vcf_valid.pl -N {wildcards.sample} -E -f {params.AF_thr} > {output} 2> {log}
#        """

#rule freebayes:
#    input:
#        bam=PATH_BAM + '/{sample}.bam'
#    output:
#        PATH_VAR+'/{sample}.vcf'
#    log:
#        PATH_LOG+ '/{sample}_freebayes'
#    params:
#        ref = config['vardict']['reference'],
#        target = config['vardict']['target']
#    conda:
#        '../env/freebayes.yaml'
#    threads:
#        config['Ncores']
#    shell:
#        """
#        freebayes-parallel <(fasta_generate_regions.py {params.ref}.fai 100000) {threads} -f {params.ref} {input.bam} --target {params.target} > {output} 2> {log}
#        """

################# Pipeline from Jinhyuk Bhin ##########################
# split target regions for freebayes parallele running
rule sort_targets:
    input:
        config["freebayes"]["targets"]
    output:
        temp(path.join(PATH_VAR, 'target', 'targets.sorted.bed'))
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


rule merge_targets:
    input:
        path.join(PATH_VAR, 'target', 'targets.sorted.bed')
    output:
        temp(path.join(PATH_VAR, 'target', 'targets.merged.bed'))
    shell:
        "bedtools merge -i {input} > {output}"

rule split_targets:
    input:
        path.join(PATH_VAR, 'target', 'targets.merged.bed')
    output:
        expand(path.join(PATH_VAR, 'target', 'targets.{part}.bed'),
        part=range(config["freebayes"]["njobs"]))
    run:
        targets = pd.read_csv(input[0], sep="\t", header=None)
        if len(targets) < len(output):
            raise ValueError("Number of freebayes threads must be less than or "
                        "equal to the number of (merged) target regions")

        chunk_size = len(targets) / len(output)
        chunks = (np.arange(len(targets)) / chunk_size).astype(int)
        for i, chunk in targets.groupby(chunks):
            chunk.to_csv(output[i], sep="\t", index=False, header=None)

rule freebayes:
    input:
        bams=expand(path.join(PATH_BAM, '{sample}.rg.bam'), sample=IDs),
        index=expand(path.join(PATH_BAM, '{sample}.rg.bam.bai'), sample=IDs),
        targets=path.join(PATH_VAR, 'target', 'targets.{part}.bed')
    output:
        temp(path.join(PATH_VAR, 'split', 'call.{part}.vcf.gz'))
    params:
        reference=config["freebayes"]["reference"],
        extra=config["freebayes"]["extra"]
    conda:
        '../env/freebayes.yaml'
    log:
        path.join(PATH_LOG, 'freebayes', 'freebayes.log')
    shell:
        "freebayes {params.extra} --fasta-reference {params.reference} "
        "--targets {input.targets} {input.bams} 2> {log} | "
        "bcftools view --output-file {output} --output-type z && "
        "bcftools index {output}"


rule combine_vcfs:
    input:
        expand(path.join(PATH_VAR, 'split', 'call.{part}.vcf.gz'),
        part=range(config["freebayes"]["njobs"]))
    conda:
        '../env/freebayes.yaml'
    output:
        path.join(PATH_VAR, 'combined.vcf.gz')
    shell:
        "bcftools concat "
        "--allow-overlaps --remove-duplicates "
        "--output-type z --output {output} {input} && "
        "bcftools index {output}"


##############################################################################


#rule annotation_vcf:
#    input:
#        path.join(PATH_VAR, 'combined.vcf.gz')
#    output:
#        path.join(PATH_VAR, 'combined.ann.vcf.gz')
#    conda:
#        '../env/vcftools.yaml'
#    log:
#        path.join(PATH_LOG, 'SnpSift_annotate.log')
#    params:
#        known=
#    shell:
#        """
#        export _JAVA_OPTIONS=-Xmx16G; SnpSift annotate {parmas.known} {input} > {output} 2> {log}
#        """

#rule filter_known:
#    input:
#        path.join(PATH_VAR, 'combined.ann.vcf.gz')
#    output:
#        path,join(PATH_VAR, 'combined.ann.filtered.vcf')
#    params:
#        "!(exists KGPhase1) && !(exists

