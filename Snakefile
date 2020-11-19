#!/usr/bin/env python3

# import peppy
# proj = peppy.Project('config/config.yaml')
# proj.get_sample('n5d2')

from pathlib import Path

#############
# FUNCTIONS #
#############

def get_bc(wildcards):
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return my_pep['bc']


def get_sample_reads(wildcards):
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    if my_pep['group'] == 'castes':
        my_key = f'l1r{wildcards.r}'
        return(my_pep[my_key])
    if my_pep['group'] == 'angry':
        my_keys = [f'l1r{wildcards.r}', f'l2r{wildcards.r}']
        return [my_pep[x] for x in my_keys]


def pick_trinity_input(wildcards):
    if wildcards.run == 'merged':
        return {
            'r1':
            expand('output/000_tmp/{sample}.r1_joined_with_merged.fastq',
                   sample=all_samples),
            'r2':
            expand('output/020_merged/{sample}_R2.fastq',
                   sample=all_samples)}
    elif wildcards.run == 'raw':
        return {
            'r1':
            expand('output/010_reads/{sample}_R1.fastq',
                   sample=all_samples),
            'r2':
            expand('output/010_reads/{sample}_R2.fastq',
                   sample=all_samples)}
    else:
        raise ValueError(f'wtf run {wildcards.run}')


def posix_path(x):
    return(Path(x).resolve().as_posix())


###########
# GLOBALS #
###########

bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
bioconductor = 'shub://TomHarrop/r-containers:bioconductor_3.11'
pandas_container = 'shub://TomHarrop/py-containers:pandas_0.25.3'
trinity = 'shub://TomHarrop/assemblers:trinity_2.11.0'
star = 'shub://TomHarrop/align-utils:star_2.7.6a'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'

# samples
pepfile: 'config/config.yaml'
all_samples = sorted(set(pep.sample_table['sample_name']))

ref = 'data/ref/GCA_014466185.1_ASM1446618v1_genomic.fna'
gff = 'data/ref/GCA_014466185.1_ASM1446618v1_genomic.gff'


#########
# RULES #
#########

rule target:
    input:
        expand('output/030_trinity/trinity.{run}/read_partitions.tar',
               run=['merged', 'raw']),
        expand('output/040_trinity-abundance/{run}/salmon.isoform.counts.matrix',
               run=['merged', 'raw']),
        expand('output/050_deseq/{run}/dds.salmon.Rds',
               run=['merged', 'raw'])


# deseq2
rule generate_deseq_object_salmon:
    input:
        quant_files = expand(
            'output/040_trinity-abundance/{{run}}/{sample}/quant.sf',
            sample=all_samples),
        tx2gene = 'output/030_trinity/trinity.{run}/Trinity.fasta.gene_trans_map'
    output:
        dds = 'output/050_deseq/{run}/dds.salmon.Rds'
    log:
        'output/logs/generate_deseq_object.{run}.log'
    singularity:
        bioconductor
    script:
        'src/generate_deseq_object.R'

# re-map reads
rule abundance_to_matrix:
    input:
        qf = expand('output/040_trinity-abundance/{{run}}/{sample}/quant.sf',
                    sample=all_samples),
        gtm = 'output/030_trinity/trinity.{run}/Trinity.fasta.gene_trans_map',
    output:
        'output/040_trinity-abundance/{run}/salmon.isoform.counts.matrix',
        'output/040_trinity-abundance/{run}/salmon.isoform.TPM.not_cross_norm'
    params:
        outdir = 'output/040_trinity-abundance/{run}',
        qf = lambda wildcards, input: ' '.join([posix_path(x) for x in input.qf])
    log:
        Path('output/logs/abundance_to_matrix.{run}.log').resolve()
    singularity:
        trinity
    shell:
        'cd {params.outdir} || exit 1 ; '
        'abundance_estimates_to_matrix.pl '
        '--est_method salmon '
        '--gene_trans_map ' + posix_path('{input.gtm}') + ' '
        '--name_sample_by_basedir '
        '--basedir_index -2 '
        '{params.qf} '
        # posix_path('{input.qf}') +
        '&> {log}'

rule trinity_abundance:
    input:
        'output/030_trinity/trinity.{run}/Trinity.fasta.salmon.idx',
        transcripts = 'output/030_trinity/trinity.{run}/Trinity.fasta',
        r1 = expand('output/010_reads/{sample}_R1.fastq',
                    sample=all_samples),
        r2 = expand('output/010_reads/{sample}_R2.fastq',
                    sample=all_samples),
        samples_txt = 'output/040_trinity-abundance/trinity_samples.txt'
    output:
        expand('output/040_trinity-abundance/{{run}}/{sample}/quant.sf',
               sample=all_samples)
    params:
        outdir = 'output/040_trinity-abundance/{run}',
    log:
        Path('output/logs/trinity_abundance.{run}.log').resolve()
    threads:
        workflow.cores
    singularity:
        trinity
    shell:
        'cd {params.outdir} || exit 1 ; '
        'align_and_estimate_abundance.pl '
        '--transcripts ' + posix_path('{input.transcripts}') + ' '
        '--seqType fq '
        '--samples_file ' + posix_path('{input.samples_txt}') + ' '
        '--est_method salmon '
        '--SS_lib_type RF '
        '--thread_count {threads} '
        '--trinity_mode '
        '&> {log}'

rule generate_samples_txt:
    input:
        samples_txt = 'config/trinity_samples.txt'
    output:
        samples_txt = 'output/040_trinity-abundance/trinity_samples.txt'
    container:
        pandas_container
    script:
        'src/generate_samples_txt.py'

rule trinity_abundance_prep:
    input:
        transcripts = 'output/030_trinity/trinity.{run}/Trinity.fasta'
    output:
        directory('output/030_trinity/trinity.{run}/Trinity.fasta.salmon.idx')
    log:
        'output/logs/trinity_abundance_prep.{run}.log'
    singularity:
        trinity
    shell:
        'align_and_estimate_abundance.pl '
        '--transcripts {input.transcripts} '
        '--est_method salmon '
        '--trinity_mode '
        '--prep_reference '
        '&> {log}'


# trinity
rule trinity_cleanup:
    input:
        'output/030_trinity/trinity.{run}/read_partitions'
    output:
        'output/030_trinity/trinity.{run}/read_partitions.tar'
    log:
        'output/logs/trinity_cleanup.{run}.log'
    shell:
        'tar -cvf '
        '{output} '
        '{input} '
        '&> {log}'

rule trinity:
    input:
        unpack(pick_trinity_input)
    output:
        'output/030_trinity/trinity.{run}/Trinity.fasta',
        'output/030_trinity/trinity.{run}/Trinity.fasta.gene_trans_map',
        temp(directory(
            'output/030_trinity/trinity.{run}/read_partitions'))
    params:
        outdir = 'output/030_trinity/trinity.{run}',
        r1 = lambda wildcards, input:
            ','.join(input.r1),
        r2 = lambda wildcards, input:
            ','.join(input.r2)
    log:
        'output/logs/trinity.{run}.log'
    threads:
        workflow.cores
    singularity:
        trinity
    shell:
        'Trinity '
        # '--FORCE '
        '--seqType fq '
        '--max_memory 800G '
        '--left {params.r1} '
        '--right {params.r2} '
        '--SS_lib_type RF '
        '--CPU {threads} '
        '--output {params.outdir} '
        '&> {log}'



# merge the input reads, try with and without
rule join_merged_with_r1:
    input:
        r1 = 'output/020_merged/{sample}_R1.fastq',
        merged = 'output/020_merged/{sample}_merged.fastq'
    output:
        temp('output/000_tmp/{sample}.r1_joined_with_merged.fastq')
    singularity:
        bbduk
    shell:
        'cat {input.r1} {input.merged} > {output}'

rule merge:
    input:
        r1 = 'output/010_reads/{sample}_R1.fastq',
        r2 = 'output/010_reads/{sample}_R2.fastq'
    output:
        merged = 'output/020_merged/{sample}_merged.fastq',
        r1 = 'output/020_merged/{sample}_R1.fastq',
        r2 = 'output/020_merged/{sample}_R2.fastq',
        ihist = 'output/020_merged/{sample}.ihist.txt'
    params:
        adaptors = '/adapters.fa'
    log:
        'output/logs/merge.{sample}.log'
    threads:
        workflow.cores
    singularity:
        bbduk
    shell:
        'bbmerge.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.merged} '
        'outu1={output.r1} '
        'outu2={output.r2} '
        'ihist={output.ihist} '
        'verystrict=t '
        'adapters={params.adaptors} '
        '2>{log}'

# # trinity doesn't do interleaved
# rule dont_split:
#     input:
#         'output/000_tmp/{indiv}.trim.fastq'
#     output:
#         r1 = 'output/010_reads/{indiv}_R1.fastq.gz'
#     wildcard_constraints:
#         indiv = '|'.join(indiv_r1_only),
#     log:
#         'output/logs/split.{indiv}.log'
#     singularity:
#         bbduk
#     shell:
#         'reformat.sh '
#         'in={input} '
#         'int=f '
#         'out={output.r1} '
#         'zl=9 '
#         '2> {log}'

rule split:
    input:
        'output/000_tmp/{sample}.trim.fastq'
    output:
        r1 = 'output/010_reads/{sample}_R1.fastq',
        r2 = 'output/010_reads/{sample}_R2.fastq'
    log:
        'output/logs/split.{sample}.log'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'reformat.sh '
        'in={input} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        '2> {log}'

rule trim:
    input:
        'output/000_tmp/{sample}.decon.fastq'
    output:
        temp('output/000_tmp/{sample}.trim.fastq')
    params:
        trim = '/adapters.fa'
    log:
        log = 'output/logs/trim.{sample}.log',
        stats = 'output/logs/trim.{sample}.stats'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '

rule decon:
    input:
        'output/000_tmp/{sample}.repair.fastq'
    output:
        pipe('output/000_tmp/{sample}.decon.fastq')
    params:
        filter = '/phix174_ill.ref.fa.gz'
    log:
        log = 'output/logs/decon.{sample}.log',
        stats = 'output/logs/decon.{sample}.stats'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '

rule repair:
    input:
        'output/000_tmp/{sample}.barcode.fastq'
    output:
        # this should be a pipe but it won't let me
        temp('output/000_tmp/{sample}.repair.fastq')
    log:
        'output/logs/repair.{sample}.log'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'repair.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'repair=t '
        '>> {output} '
        '2> {log}'

rule check_barcodes:
    input:
        r1 = 'output/000_tmp/{sample}.r1.fastq',
        r2 = 'output/000_tmp/{sample}.r2.fastq'
    output:
        pipe('output/000_tmp/{sample}.barcode.fastq')
    params:
        bc = lambda wildcards: get_bc(wildcards)
    log:
        'output/logs/check_barcodes.{sample}.log'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'reformat.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'int=f '
        'out=stdout.fastq '
        'barcodefilter=t '
        'barcodes={params.bc} '
        '>> {output} '
        '2> {log}'

# make this temp later, but we'll use it for fastqc
rule combine_read_file:
    input:
        get_sample_reads
    output:
        'output/000_tmp/{sample}.r{r}.fastq'
    shell:
        'zcat {input} > {output}'


# STAR mapping for genome-guided mode
rule star_target:
    input:
        'output/030_trinity/trinity.genome/read_partitions.tar'


rule trinity_genome:
    input:
        bam = 'output/025_star/merged.bam'
    output:
        'output/030_trinity/trinity.genome/Trinity.fasta',
        'output/030_trinity/trinity.genome/Trinity.fasta.gene_trans_map',
        temp(directory(
            'output/030_trinity/trinity.genome/read_partitions'))
    params:
        outdir = 'output/030_trinity/trinity.genome',
    log:
        'output/logs/trinity.genome.log'
    threads:
        workflow.cores
    singularity:
        trinity
    shell:
        'Trinity '
        # '--FORCE '
        '--max_memory 800G '
        ' --genome_guided_bam {input.bam} '
        ' --genome_guided_max_intron 10000 '
        '--SS_lib_type RF '
        '--CPU {threads} '
        '--output {params.outdir} '
        '&> {log}'


rule merge_bam:
    input:
        bam = expand('output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
                     sample=all_samples)
    output:
        'output/025_star/merged.bam'
    log:
        'output/logs/merge_bam.log'
    threads:
        min(20, workflow.cores)
    container:
        samtools
    shell:
        'samtools merge '
        '-l 9 '
        '-O BAM '
        '-@ {threads} '
        '{output} '
        '{input.bam} '
        '2> {log}'


rule star_second_pass:
    input:
        r1 = 'output/010_reads/{sample}_R1.fastq',
        r2 = 'output/010_reads/{sample}_R2.fastq',
        star_reference = 'output/007_star-index/SA',
        junctions = expand('output/025_star/pass1/{sample}.SJ.out.tab',
                           sample=all_samples)
    output:
        bam = 'output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
        counts = 'output/025_star/pass2/{sample}.ReadsPerGene.out.tab'
    threads:
        workflow.cores
    params:
        genome_dir = 'output/007_star-index',
        prefix = 'output/025_star/pass2/{sample}.'
    log:
        'output/logs/star_second_pass.{sample}.log'
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outReadsUnmapped Fastx '
        '--quantMode GeneCounts '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_first_pass:
    input:
        r1 = 'output/010_reads/{sample}_R1.fastq',
        r2 = 'output/010_reads/{sample}_R2.fastq',
        star_reference = 'output/007_star-index/SA'
    output:
        sjdb = 'output/025_star/pass1/{sample}.SJ.out.tab'
    threads:
        workflow.cores
    params:
        genome_dir = 'output/007_star-index',
        prefix = 'output/025_star/pass1/{sample}.'
    log:
        'output/logs/star_first_pass.{sample}.log'
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '          # troubleshoot gtf
        # '--outSAMtype SAM '               # troubleshoot gtf
        # '--quantMode GeneCounts '       # troubleshoot gtf
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_index:
    input:
        fasta = ref,
        gff = gff
    output:
        'output/007_star-index/SA'
    params:
        outdir = 'output/007_star-index'
    log:
        'output/logs/star_index.log'
    threads:
        workflow.cores
    container:
        star
    shell:
        'STAR '
        'runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gff} '
        '--genomeSAindexNbases 12 '
        '--sjdbGTFtagExonParentTranscript Parent '
        '--sjdbGTFtagExonParentGene locus_tag '
        '&> {log}'


