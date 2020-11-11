#!/usr/bin/env python3

from pathlib import Path
import pandas

# sample_table_file = 'config/trinity_samples.txt'

sample_table_file = snakemake.input['samples_txt']

sample_table = pandas.read_csv(
    sample_table_file,
    sep='\t',
    header=None)

sample_table[2] = [Path(x).resolve().as_posix() for x in sample_table[2]]
sample_table[3] = [Path(x).resolve().as_posix() for x in sample_table[3]]

sample_table.to_csv(
    snakemake.output['samples_txt'],
    sep='\t',
    header=False,
    index=False)
