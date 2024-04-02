from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from tqdm import tqdm
import numpy as np
import pickle
import argparse
import pandas as pd
from collections import Counter
from itertools import islice
import os

parser = argparse.ArgumentParser(description='Retrieve all barcodes found in pair-end fastq files')

parser.add_argument('--sample', type=str, default='',
                    help='name of the sample')

parser.add_argument('--suffix', type=str, default='pear',
                    help='suffix of output')

parser.add_argument('--datadir', type=str, default='./fastq',
                    help='location of fastq files')

parser.add_argument('--savedir', type=str, default='',
                    help='path to save dir')

parser.add_argument('--pear_output', type=str, default=".",
                    help='path to save PEAR-merged fastq output')

parser.add_argument('--minphred', type=int, default=25,
                    help='min phred score')

parser.add_argument('--offset_primer', type=int, default=5,
                    help='number of offset nucleotides for mapping')

parser.add_argument('--PCR_primer_F', type=str, default="GGATTATTCATACCGTCCCA",
                    help='forward PCR primer')

parser.add_argument('--PCR_primer_R', type=str, default="CAAATGTGGTATGGCTGATT",
                    help='reverse PCR primer')

args = parser.parse_args()
sample = args.sample
suffix = args.suffix
datadir = args.datadir
pear_output = args.pear_output
min_phred = args.minphred
savedir = args.savedir
offset_primer = args.offset_primer
PCR_primer_F = args.PCR_primer_F
PCR_primer_R = args.PCR_primer_R

fastbac = [PCR_primer_F, str(Seq(PCR_primer_R).reverse_complement())]
fastbac_r = [PCR_primer_R, str(Seq(PCR_primer_F).reverse_complement())]

def extract_bc(seq, primer, offset=offset_primer):
    """
    extract barcode sequence using PCR primer sequences
    """
    primer = [primer[0][offset:], primer[1][:-offset]]
    if (primer[0] in seq) and (primer[1] in seq):
        start = seq.index(primer[0]) + len(primer[0])
        end = seq.index(primer[1])
        return seq[start:end], start, end
    else:
        return None
        
def merge_read_rev(s, primers = fastbac_r):
    '''
    inputs: PEAR-merged read
        reverse read will have (fastbac_r)
    '''
    if extract_bc(s.seq, primers):
        bc, st, e = extract_bc(s.seq, primers)
        phred = s.letter_annotations['phred_quality'][st:e]
        if (len(bc)>0) and (len(phred)>0):
            return str(Seq(bc).reverse_complement()), phred
        else:
            return None
    else:
        return None
    
def merge_read_for(s, primers = fastbac):
    '''
    inputs: PEAR-merged read
        forward read will have (fastbac)
    '''
    if extract_bc(s.seq, primers):
        bc, st, e = extract_bc(s.seq, primers)
        phred = s.letter_annotations['phred_quality'][st:e]
            
        if (len(bc)>0) and (len(phred)>0):
            return str(bc), phred
        else:
            return None
    else:
        return None

def retrieve_barcode(sample, phred_threshold, datadir):

    fastq = SeqIO.parse(f"{datadir}/{sample}.assembled.fastq", 'fastq')

    barcodes = []
    total_reads = 0                     # total PE reads after PEAR
    mapped_reads = 0                    # reads containing PCR primers
    mapped_highquality_reads = 0        # reads containing PCR primers and high Phred score
    unmapped_reads = 0                  # unidentified reads after merge 2 PE reads

    for s in tqdm(fastq):
        total_reads += 1

        # if forward reads
        if merge_read_for(s):
            mapped_reads += 1
            try:
                bc, phred = merge_read_for(s)
                if min(phred) >= phred_threshold:
                    barcodes.append(bc)
                    mapped_highquality_reads += 1
            except Exception:
                print(merge_read_for(s))

        # if reverse reads
        elif merge_read_rev(s):
            mapped_reads += 1
            try:
                bc, phred = merge_read_rev(s)
                if min(phred) >= phred_threshold:
                    barcodes.append(bc)
                    mapped_highquality_reads += 1
            except Exception:
                print(merge_read_rev(s))

        else:
            unmapped_reads += 1

    assert total_reads == mapped_reads + unmapped_reads
    summary_count = {
        "total_reads": [total_reads],
        "mapped_reads": [mapped_reads],
        "mapped_highquality_reads": [mapped_highquality_reads],
        "unmapped_reads": [unmapped_reads],
    }
    return barcodes, summary_count
    
####################################################################################
####################################################################################
####################################################################################

# create subfolder to save PEAR output
os.system(f"mkdir {pear_output}/{sample}")

# merge paried-end reads
cmd = f"pear \
-f {datadir}/{sample}_R1_001.fastq.gz \
-r {datadir}/{sample}_R2_001.fastq.gz \
-o {pear_output}/{sample}/{sample} \
-j 32 \
-b {min_phred} -q 15 > {savedir}/PEAR_{sample}.log"

os.system(cmd)

# count barcode from merged fastq files
barcodes, summary_count = retrieve_barcode(sample, phred_threshold = min_phred, datadir=f"{pear_output}/{sample}")
barcodes_count = Counter(barcodes)

# remove the generated fastq files to save space
os.system(f"rm -r {pear_output}/{sample}")

# save files
df = pd.DataFrame(list(barcodes_count.items()), columns=['barcode', 'n_bc'])
df.to_csv(f'{savedir}/barcodes_{sample}_{suffix}.csv', index=False)

df_sum = pd.DataFrame(summary_count)
df_sum.to_csv(f'{savedir}/log_{sample}_{suffix}.csv', index=False)