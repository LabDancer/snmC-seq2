#!/bin/bash 

module load BEDTools/2.30.0-GCC-11.2.0 

for i in /data/gent/vo/000/gvo00056/vsc42347/InternshipLeuven/H1_fastq/11_binning/4samples/*.gz
do
    s=$(basename $f .bam)
    bedtools coverage -a hg19.24.10kb.part -b $i > ./"$s"_coverage_10kb.bed
    
done

