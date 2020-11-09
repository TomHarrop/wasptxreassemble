#!/usr/bin/env python3

# import peppy
# proj = peppy.Project('config/config.yaml')
# proj.get_sample('n5d2')


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


###########
# GLOBALS #
###########

bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
trinity = 'shub://TomHarrop/assemblers:trinity_2.11.0'


# samples
pepfile: 'config/config.yaml'
all_samples = sorted(set(pep.sample_table['sample_name']))


#########
# RULES #
#########

rule target:
    input:
        expand('output/030_trinity/trinity.{run}/read_partitions.tar',
               run=['merged', 'raw'])

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

