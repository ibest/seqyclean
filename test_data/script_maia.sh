#!/bin/bash
FILES=/home/maia/CysticFibrosis/666CCAAXX_rRNAdepleted/Sample_HAP94_Flagellin/*.fastq.gz

for f in $FILES
do
  seqyclean -qual 24 24 -w0 20 -w1 5 -minimum_read_length 50 --new2old_illumina -U $f -o cleaned_w0.20.w1.5_50
done
