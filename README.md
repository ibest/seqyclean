# Description

Program ```SeqyClean```
Version: ```1.9.7 (2014-11-01)```

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
usage: ./seqyclean libflag input_file_name_1 [libflag input_file_name_2] -o output_prefix [options]
```

The parameter ```libflag``` here is a library type: -454 for Roche 454 reads, -1, -2 for paired-end Illumina reads, -U for single-end reads. See examples below.
            
## Common arguments for all library types
```
-h, --help - Show this help and exit.
-v <filename> - Turns on vector trimming, default=off. <filename> - is a path to a FASTA-file containing vector genomes.
-c <filename> - Turns on contaminants screening, default=off, <filename> - is a path to a FASTA-file containing contaminant genomes.
-k <value> - Common size of k-mer, default=15.
-d - Distance between consecutive k-mers, default=1.
-kc <value> - Size of k-mer used in sampling contaminat genome, default=15.
-qual <value> <value> - Turns on quality trimming, default=off. Error boundaries: max_average_error (default=20), max_error_at_ends (default=20).
-qual_only - Performs only quality trimming without trimming of adapters, default=off.
-ow - Overwrite existing results, default=off.
-minlen <value> - Minimum length of read to accept, default=50 bp.
-polyat [cdna] [cerr] [crng] - Turns on poly A/T trimming, default=off. Parameters: cdna (default=10) - maximum size of a poly tail, cerr (default=3) - maximum number of G/C nucleotides within a tail, cnrg (default=50) - range to look for a tail within a read.
-verbose - Verbose output, default=off.
-detrep - Generate detailed report for each read, default=off.
```
## Roche 454 arguments
```
-t <value> - Number of threads (not yet applicable to Illumina mode), default=4.
-fastq - Output in FASTQ format, default=off.
-fasta_out - Output in FASTA format, default=off.
-m <filename> - Using custom barcodes, default=off. <filename> - a path to a FASTA-file with custom barcodes.
```
## Illumina paired- and single-end arguments
```
-1 <filename1> -2 <filename2> - Paired-end mode (see examples below)
-U <filename> - Single-end mode
-shuffle - Store non-paired Illumina reads in shuffled file, default=off.
-i64 - Turns on 64-quality base, default = off.
-adp <filename> - Turns on using custom adapters, default=off. <filename> - FASTA file with adapters
-alen <value> - Minimum adapter length for dovetail overlap for adapter trimming, default = 60 bp.
-at <value> - overlap threshold (only in paired-end mode, default = 0.75.
-overlap <value> - Turns on merging overlapping paired-end reads. Parameter <value> is a minimum overlap length, default=16 bp.
-dup [-startdw][-sizedw] - Turns on screening duplicated sequences, default=off. Here startdw (defalt=10) and sizedw (default=15) are starting position and size of the window within a read.
-no_ts_adapter_trim - Turn off TruSeq adapters trimming, default=off.
-new2old - A switch to fix read IDs, default=off ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 ).
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

# Contacts

Please ask Ilya (ilyaz@uidaho.edu or zhba3458@vandals.uidaho.edu) in case of any questions.