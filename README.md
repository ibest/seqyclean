# Description

Program ```SeqyClean```
Version: ```1.9.6 (2014-10-28)```

Main purpose of this software is to pre-process NGS data in order to prepare for downstream analysis.

SeqyClean offers:

* Adapter/key/primers filtering
* Vector and contaminants filtering.
* Quality trimming.
* Poly A/T trimming.
* Overlapping paired reads.

It handles SFF and FASTQ file formats.

# Usage  
```
usage: ./seqyclean libflag input_file_name_1 [libflag input_file_name_2] -o output_prefix [options]\n"
```            
## Common arguments for all library types
```
-h, --help - show this help and exit.\n"
-v <filename> - Turns on vector trimming, default=off. <filename> - is a path to a FASTA-file containing vector genomes.\n"
-c <filename> - Turns on contaminants screening, default=off, <filename> - is a path to a FASTA-file containing contaminant genomes.\n"
-k <int> - Common size of k-mer, default=15\n"
-d - distance between consecutive k-mers, default=1\n"
-kc <int> - Size of k-mer used in sampling contaminat genome, default=15\n"
-qual <int> <int> - Turns on quality trimming, default=off. Error boundaries: max_average_error (default=20), max_error_at_ends (default=20)\n"
-qual_only - Performs only quality trimming without trimming of adapters, default=off.\n"
-ow - Overwrite existing results, default=off\n"
-minlen <int> - Minimum length of read to accept, default=50 bp.\n"
-polyat [cdna] [cerr] [crng] - Turns on poly A/T trimming, default=off. Parameters: cdna (default=10) - maximum size of a poly tail, cerr (default=3) - maximum number of G/C nucleotides within a tail, cnrg (default=50) - range to look for a tail within a read.\n"
-verbose - Verbose output, default=off.\n"
-detrep - Generate detailed report for each read, default=off.\n"
```
## Roche 454 arguments
```
-t <int> - Number of threads (not yet applicable to Illumina mode), default=4.\n" 
-fastq - Output in FASTQ format, default=off.\n"
-fasta_out - Output in FASTA format, default=off.\n"
-m <filename> - Using custom barcodes, default=off. <filename> - a path to a FASTA-file with custom barcodes.\n"
```
## Illumina paired- and single-end arguments
```
-1 <filename1> -2 <filename2> - Paired-end mode (see examples below)
-U <filename> - single-end mode
-shuffle - Store non-paired Illumina reads in shuffled file, default=off.\n"
-i64 - Turns on 64-quality base, default = off.\n"
-adp <filename> - Turns on using custom adapters, default=off. <filename> - FASTA file with adapters\n"
-alen <value> - minimum adapter length for dovetail overlap, default = 60 bp.\n"
-ot <value> - overlap threshold (only in paired-end mode, default = 0.75.\n"
-dup [-startdw][-sizedw] - Turns on screening duplicated sequences, default=off. Here startdw (defalt=10) and sizedw (default=15) are starting position and size of the window within a read.\n" 
-no_ts_adapter_trim - Turn off TruSeq adapters trimming, default=off.\n"
-new2old - switch to fix read IDs, default=off ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 ).\n";
```

# Examples

## Roche 454
Output in SFF, no quality trimming, vector trimming is performed:
```
./seqyclean -454 test_data/in.sff -o test/Test454 -v test_data/vectors.fasta
```
Output in SFF, quality trimming with default parameters, vector trimming and contaminants screening are performed:
```
./seqyclean -454 test_data/in.sff -o test/Test454 -qual -v test_data/vectors.fasta -c test_data/contaminants.fasta
```

## Illumina

### Paired-end
Trimming of adapters is performed, quality trimming with default parameters:
```
./seqyclean -1 test_data/R1.fastq.gz -2 test_data/R2.fastq.gz -qual -o test/Test_Illumina
``` 
   
Trimmings of adapters and vectors are performed, quality trimming with default parameters:
```
./seqyclean -1 test_data/R1.fastq.gz -2 test_data/R2.fastq.gz -qual -v test_data/vectors.fasta -o test/Test_Illumina
```    

### Single-end
Trimming of adapters, vectors and contaminant screening are performed, quality trimming with default parameters:
```
./seqyclean -U test_data/R1.fastq.gz -o test/Test_Illumina -v test_data/vectors.fasta -c test_data/contaminants.fasta
                                
```


# Supported RL MIDS (for Roche 454 only)

```
RL1,ACACGACGACT,RL1,AGTCGTGGTGT
RL2,ACACGTAGTAT,RL2,ATACTAGGTGT
RL3,ACACTACTCGT,RL3,ACGAGTGGTGT
RL4,ACGACACGTAT,RL4,ATACGTGGCGT
RL5,ACGAGTAGACT,RL5,AGTCTACGCGT
RL6,ACGCGTCTAGT,RL6,ACTAGAGGCGT
RL7,ACGTACACACT,RL7,AGTGTGTGCGT
RL8,ACGTACTGTGT,RL8,ACACAGTGCGT
RL9,ACGTAGATCGT,RL9,ACGATCTGCGT
RL10,ACTACGTCTCT,RL10,AGAGACGGAGT
RL11,ACTATACGAGT,RL11,ACTCGTAGAGT
RL12,ACTCGCGTCGT,RL12,ACGACGGGAGT
RL13,AGACTCGACGT,RL13,ACGTCGGGTCT
RL14,AGTACGAGAGT,RL14,ACTCTCGGACT
RL15,AGTACTACTAT,RL15,ATAGTAGGACT
RL16,AGTAGACGTCT,RL16,AGACGTCGACT
RL17,AGTCGTACACT,RL17,AGTGTAGGACT
RL18,AGTGTAGTAGT,RL18,ACTACTAGACT
RL19,ATAGTATACGT,RL19,ACGTATAGTAT
RL20,CAGTACGTACT,RL20,AGTACGTGCTG
RL21,CGACGACGCGT,RL21,ACGCGTGGTCG
RL22,CGACGAGTACT,RL22,AGTACTGGTCG
RL23,CGATACTACGT,RL23,ACGTAGTGTCG
RL24,CGTACGTCGAT,RL24,ATCGACGGACG
RL25,CTACTCGTAGT,RL25,ACTACGGGTAG
RL26,GTACAGTACGT,RL26,ACGTACGGTAC
RL27,GTCGTACGTAT,RL27,ATACGTAGGAC
RL28,GTGTACGACGT,RL28,ACGTCGTGCAC
RL29,ACACAGTGAGT,RL29,ACTCACGGTGT
RL30,ACACTCATACT,RL30,AGTATGGGTGT
RL31,ACAGACAGCGT,RL31,ACGCTGTGTGT
RL32,ACAGACTATAT,RL32,ATATAGTGTGT
RL33,ACAGAGACTCT,RL33,AGAGTCTGTGT
RL34,ACAGCTCGTGT,RL34,ACACGAGGTGT
RL35,ACAGTGTCGAT,RL35,ATCGACAGTGT
RL36,ACGAGCGCGCT,RL36,AGCGCGCGCGT
```

# Contacts

Please ask Ilya (zhba3458@vandals.uidaho.edu) in case of any questions.