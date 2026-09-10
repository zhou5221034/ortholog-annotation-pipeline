#!/bin/bash
if [ $# -lt 2 ]; then
    echo "Usage：sh Getpep.sh [genome.fa] [gff]"
    exit 1
fi
genome=$1
gff=$2

perl LongestGff.pl $gff
gffread -g $genome -x out.cds $gff
perl cds2aa.pl out.cds > out.pep
perl GetSeq.pl mapid out.pep > out.longest.pep

echo `date`
