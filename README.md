# Count barcodes from amplicon paired-end sequencing data

retrieve_barcode.py script will run directly from paired-end fastq files (with suffix '_R1_001.fastq.gz' and '_R2_001.fastq.gz'). \n
retrieve_barcode_PEAR.py script will instead merge 2 paired-end files into one by PEAR, and then count. It's much faster and gives generally higher barcode numbers.
