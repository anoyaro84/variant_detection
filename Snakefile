

import pandas as pd
import os
import glob
from os import path as path
import numpy as np

configfile: 'config.yaml'

# take information from config
PATH_FASTQ = config["path"]["fastq"]
PATH_VAR = config["path"]["variant"]
PATH_QC = config["path"]["qc"]
PATH_LOG = config["path"]["log"]
PATH_BAM = config['path']['bam']
START_BAM = config['Use_BAM']

# some options are taken from environmental variables
if config['Use_global'] is True:
    for env in config['environmental']:
        if config['environmental'][env] in os.environ:
            globals()[env] = os.environ[config['environmental'][env]]
            print(env + ' is overrived by enviromental variable ' + config['environmental'][env] + ': ' + globals()[env])


# if samples table exists, take fastq listed in the table,
# otherwise, it will just take every fastqs in PATH_FASTQ
if os.path.isfile('samples.tsv'):
    samples = pd.read_csv('samples.tsv', sep='\t')
    IDs = samples.IDs
    print("obtaining samples from config.yaml")
elif not START_BAM:
    print("obtaining samples from the path : " + PATH_FASTQ)
    Temp = glob.glob(PATH_FASTQ+"/*.fastq.gz")
    IDs = list(set([os.path.basename(tmp).split('_')[0] for tmp in Temp]))
else:
    print("obtaining aligned samples from the path : " + PATH_BAM)
    Temp = glob.glob(PATH_BAM+'/*.bam')
    IDs = [os.path.splitext(os.path.basename(tmp))[0] for tmp in Temp]

print('found ' + str(len(IDs)) + ' files to be processed')

rule all:
    input:
        PATH_QC+'/multiqc.html',
        PATH_VAR + "/combined.vcf.gz",

include: "rules/variant.snakefile"
