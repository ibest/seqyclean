# Description

Program ```SeqyClean```
Version: ```1.9.1 (2014-04-13)```

Main purpose of this software is to clean reads. It provide adapter/key/primers searching and quality trimming (LUCY).

# Usage  
## Roche 454

```
./seqyclean -454 input_file_name -o output_prefix [-v vector_file]
						  [-c file_of_contaminants]
						  [-m file_of_RL_MIDS]
						  [-k k_mer_size]
						  [-kc k_mer_size]
						  [-f overlap ]
						  [-t number_of_threads]
						  [-qual max_avg_error max_error_at_ends]
						  [--qual_only]
						  [--fastq]
						  [--keep_fastq_orig]
                                                  [--ow]
						  [-minimum_read_length <value>]
						  [-polyat [cdna] [cerr] [crng] ]
```

## Illumina

### Paired-end
```
./seqyclean -1 input_file_name_1 -2 input_file_name_2 -o output_prefix [-v vector_file]
                                                                       [-c file_of_contaminants]
                                                                       [-k k_mer_size]
                                                                       [-kc k_mer_size]
                                                                       [-qual max_avg_error max_error_at_ends]
                                                                       [--qual_only]
                                                                       [-minimum_read_length <value>]
                                                                       [-i64]
                                                                       [-adapter_length <value>]
                                                                       [-ot <value>]
                                                                       [--overlap <minoverlap=value>]
                                                                       [--ow]
                                                                       [--dup][-start_dw][-size_dw]
                                                                       [-polyat [cdna] [cerr] [crng] ]
                                                                       [--no_ts_adapter_trim ]
                                                                       [--new2old_illumina] - switch to fix read IDs ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 )
```

### Single-end

```
./seqyclean -U input_file_name -o output_prefix [-v vector_file]
                                                                       [-c file_of_contaminants]
                                                                       [-k k_mer_size]
                                                                       [-kc k_mer_size]
                                                                       [-qual max_avg_error max_error_at_ends]
                                                                       [--qual_only]
                                                                       [-minimum_read_length <value>]
                                                                       [--ow]
                                                                       [-polyat [cdna] [cerr] [crng] ]
                                                                       [--no_ts_adapter_trim ]
                                                                       [--new2old_illumina] - switch to fix read IDs ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 )
```


# Examples

## Roche 454
```
./seqyclean -454 test_data/in.sff -o test/Test454 -v test_data/vectors.fasta
```

## Illumina
```
./seqyclean -1 test_data/R1.fastq.gz -2 test_data/R2.fastq.gz -o test/Test_Illumina
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

# Contact

Please ask Ilya (zhba3458@vandals.uidaho.edu) in case of any questions.


