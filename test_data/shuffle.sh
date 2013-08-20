#! /bin/bash    
while true; do
  read -r lineA <&3
  read -r lineB <&4

  if [ -z "$lineA" -o -z "$lineB" ]; then
    break
  fi
  
  echo $lineA >> shuffled.fastq
  read -r lineA <&3
  echo $lineA >> shuffled.fastq
  read -r lineA <&3
  echo $lineA >> shuffled.fastq
  read -r lineA <&3
  echo $lineA >> shuffled.fastq
  
  echo $lineB >> shuffled.fastq
  read -r lineB <&4
  echo $lineB >> shuffled.fastq
  read -r lineB <&4
  echo $lineB >> shuffled.fastq
  read -r lineB <&4
  echo $lineB >> shuffled.fastq

done 3<R1.fastq 4<R2.fastq
