# Description

Program ```SeqyClean```
Version: ```1.10.06 (2017-12-13)```

Main purpose of this software is to pre-process NGS data in order to prepare for downstream analysis.

SeqyClean offers:

* Adapter/key/primers filtering
* Vector and contaminants filtering.
* Quality trimming.
* Poly A/T trimming.
* Overlapping paired reads.

It handles SFF and FASTQ file formats.

# Prerequisites

Developer version of the ```zlib```:

```
$sudo apt-get install zlib1g-dev
```

# Installation

Clone or download the repository. Then ```cd``` to seqyclean home folder, and type ```make```.


Note: by default, it builds the binary for OS-X. It should build on Linux as well. If not, try to use this command:

```make PLATFORM=-DLINUX```

or simply contact me: ilya.zhbannikov@duke.edu

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
-qual ```max_average_error max_error_at_ends``` - Turns on quality trimming, default=off. Error boundaries: max_average_error (default=20 Phred), max_error_at_ends (default=20 Phred).
-braket ```window_size max_average_error``` - Parameter for quality trimming. By default window_size=10 and max_average_error=0.794.
-window ```window_size max_average_error``` [```window_size maximum_average_error``` [...]] - Parameters for quality trimming. By default there are two windows with size of 50 and 10 bp with the same max_average_error=0.794.
-ow - Overwrite existing results, default=off.
-minlen <value> - Minimum length of read to accept, default=100 bp.
-polyat [cdna] [cerr] [crng] - Turns on poly A/T trimming, default=off. Parameters: cdna (default=10) - maximum size of a poly tail, cerr (default=3) - maximum number of G/C nucleotides within a tail, cnrg (default=50) - range to look for a tail within a read.
-verbose - Verbose output, default=off.
-detrep - Generate detailed report for each read, default=off.
-dup [-startdw][-sizedw][-maxdup] - Turns on screening duplicated sequences, default=off. Here startdw (defalt=10) and sizedw (default=35) are starting position and size of the window within a read, -maxdup (default=3) - maximum number of duplicated sequences allowed.
-no_adapter_trim - Turns off adapter trimming, default=off.
```

## Roche 454 arguments
```
-t <value> - Number of threads (not yet applicable to Illumina mode), default=4.
-fastq - Output in FASTQ format, default=off.
-fasta - Output in FASTA format, default=off.
-m <filename> - Using custom barcodes, default=off. <filename> - a path to a FASTA-file with custom barcodes.
```
## Illumina paired- and single-end arguments
```
-1 <filename1> -2 <filename2> - Paired-end mode (see examples below)
-U <filename> - Single-end mode
-shuffle - Store non-paired Illumina reads in shuffled file, default=off.
-i64 - Turns on 64-quality base, default = off.
-adp <filename> - Turns on using custom adapters, default=off. <filename> - FASTA file with adapters
-at <value> - This option sets the similarity threshold for adapter trimming by overlap (only in paired-end mode). By default its value is set to 0.75.
-overlap <value> - This option turns on merging overlapping paired-end reads and <value> is the minimum overlap length. By default the minimum overlap length is 16 base pairs.
-new2old - A switch to fix read IDs, default=off ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 ).
-gz - A flag that indicates compressed (.gz) output, default=off.
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

# Citing SeqyClean

## BibText
```
@inproceedings{Zhbannikov:2017:SPH:3107411.3107446,
 author = {Zhbannikov, Ilya Y. and Hunter, Samuel S. and Foster, James A. and Settles, Matthew L.},
 title = {SeqyClean: A Pipeline for High-throughput Sequence Data Preprocessing},
 booktitle = {Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics},
 series = {ACM-BCB '17},
 year = {2017},
 isbn = {978-1-4503-4722-8},
 location = {Boston, Massachusetts, USA},
 pages = {407--416},
 numpages = {10},
 url = {http://doi.acm.org/10.1145/3107411.3107446},
 doi = {10.1145/3107411.3107446},
 acmid = {3107446},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {data preprocessing, high-throughput dna sequencing, sequence analysis},
} 
```

## Plain text (ACM Ref)
```
Ilya Y. Zhbannikov, Samuel S. Hunter, James A. Foster, and Matthew L. Settles. 2017. SeqyClean: A Pipeline for High-throughput Sequence Data Preprocessing. In Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics (ACM-BCB '17). ACM, New York, NY, USA, 407-416. DOI: https://doi.org/10.1145/3107411.3107446
```

## EndNote
```
%0 Conference Paper
%1 3107446
%A Ilya Y. Zhbannikov
%A Samuel S. Hunter
%A James A. Foster
%A Matthew L. Settles 
%T SeqyClean: A Pipeline for High-throughput Sequence Data Preprocessing
%B Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics
%@ 978-1-4503-4722-8
%C Boston, Massachusetts, USA
%P 407-416
%D 2017
%R 10.1145/3107411.3107446
%I ACM
```

# Contacts

Please ask Ilya (ilya.zhbannikov@duke.edu) in case of any questions.
