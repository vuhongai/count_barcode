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

parser = argparse.ArgumentParser(description='Retrieve all barcodes found in pair-end fastq files')

parser.add_argument('--sample', type=str, default='',
                    help='name of the sample')

parser.add_argument('--suffix', type=str, default='',
                    help='suffix of output')

parser.add_argument('--datadir', type=str, default='',
                    help='location of fastq files')

parser.add_argument('--savedir', type=str, default='',
                    help='path to save dir')

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
min_phred = args.minphred
savedir = args.savedir
offset_primer = args.offset_primer
PCR_primer_F = args.PCR_primer_F
PCR_primer_R = args.PCR_primer_R

fastbac = [PCR_primer_F, str(Seq(PCR_primer_R).reverse_complement())]
fastbac_r = [PCR_primer_R, str(Seq(PCR_primer_F).reverse_complement())]

def extract_bc(seq, primer, offset=offset_primer):
    primer = [primer[0][offset:], primer[1][:-offset]]
    if (primer[0] in seq) and (primer[1] in seq):
        start = seq.index(primer[0]) + len(primer[0])
        end = seq.index(primer[1])
        return seq[start:end], start, end
    else:
        return None

def merge_nt(n1,n2):
    '''
    n1=('A', 38), n2=('A', 27)
    '''
    n=0
    if n1[0]==n2[0]=='N':
        return n1
    else:
        if n1[1]>n2[1]:
            return n1
        else:
            return n2
        
def merge_read_rev(s1,s2, primers = [fastbac_r, fastbac]):
    '''
    inputs: read from fastq paired-end reads
        reverse read will have (fastbac_r in read1) and (fastbac in read2)
    '''
    if extract_bc(s1.seq, primers[0]) and extract_bc(s2.seq, primers[1]):

        bc1, str1, e1 = extract_bc(s1.seq, primers[0])
        bc2, str2, e2 = extract_bc(s2.seq, primers[1])
        
        if len(bc1)==len(bc2):
            phred1 = s1.letter_annotations['phred_quality'][str1:e1]
            phred2 = s2.letter_annotations['phred_quality'][str2:e2]

            bc2, phred2 = Seq(bc2).reverse_complement(), phred2[::-1]

            bc = [merge_nt(n1,n2)[0] for n1,n2 in zip(zip(bc1, phred1), zip(bc2, phred2))]
            phred = [merge_nt(n1,n2)[1] for n1,n2 in zip(zip(bc1, phred1), zip(bc2, phred2))]
            bc = "".join(bc)
            if (len(bc)>0) and (len(phred)>0):
                return str(Seq(bc).reverse_complement()), phred
            else:
                return None
        
        else:
            return None
    else:
        return None
    
def merge_read_for(s1,s2, primers = [fastbac, fastbac_r]):
    '''
    inputs: read from fastq paired-end reads
        reverse read will have (fastbac in read1) and (fastbac_r in read2)
    '''
    if extract_bc(s1.seq, primers[0]) and extract_bc(s2.seq, primers[1]):

        bc1, str1, e1 = extract_bc(s1.seq, primers[0])
        bc2, str2, e2 = extract_bc(s2.seq, primers[1])
        
        if len(bc1)==len(bc2):
            phred1 = s1.letter_annotations['phred_quality'][str1:e1]
            phred2 = s2.letter_annotations['phred_quality'][str2:e2]

            bc2, phred2 = Seq(bc2).reverse_complement(), phred2[::-1]

            bc = [merge_nt(n1,n2)[0] for n1,n2 in zip(zip(bc1, phred1), zip(bc2, phred2))]
            phred = [merge_nt(n1,n2)[1] for n1,n2 in zip(zip(bc1, phred1), zip(bc2, phred2))]
            bc = "".join(bc)
            
            if (len(bc)>0) and (len(phred)>0):
                return str(bc), phred
            else:
                return None
        
        else:
            return None
    else:
        return None

def retrieve_barcode(sample, phred_threshold, datadir):
    forward_reads_file = f'{datadir}/{sample}_R1_001.fastq.gz'
    reverse_reads_file = f'{datadir}/{sample}_R2_001.fastq.gz'

    with gzip.open(forward_reads_file, 'rt') as f, gzip.open(reverse_reads_file, 'rt') as r:
        forward = SeqIO.parse(f, 'fastq')
        reverse = SeqIO.parse(r, 'fastq')

        barcodes = []
        total_reads = 0
        mapped_reads = 0 # reads containing PCR primers
        mapped_highquality_reads = 0 # reads containing PCR primers and high Phred score
        unmapped_reads = 0 # unidentified reads after merge 2 PE reads

        for s1,s2 in tqdm(zip(forward,reverse)):
            total_reads += 1

            # if forward reads
            if merge_read_for(s1,s2):
                mapped_reads += 1
                try:
                    bc, phred = merge_read_for(s1,s2)
                    if min(phred) >= phred_threshold:
                        barcodes.append(bc)
                        mapped_highquality_reads += 1
                except Exception:
                    print(merge_read_for(s1,s2))

            # if reverse reads
            elif merge_read_rev(s1,s2):
                mapped_reads += 1
                try:
                    bc, phred = merge_read_rev(s1,s2)
                    if min(phred) >= phred_threshold:
                        barcodes.append(bc)
                        mapped_highquality_reads += 1
                except Exception:
                    print(merge_read_rev(s1,s2))

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
    

barcodes, summary_count = retrieve_barcode(sample, phred_threshold = min_phred, datadir=datadir)
barcodes_count = Counter(barcodes)

df = pd.DataFrame(list(barcodes_count.items()), columns=['barcode', 'n_bc'])
df.to_csv(f'{savedir}/barcodes_{sample}_{suffix}.csv', index=False)

df_sum = pd.DataFrame(summary_count)
df_sum.to_csv(f'{savedir}/log_{sample}_{suffix}.csv', index=False)