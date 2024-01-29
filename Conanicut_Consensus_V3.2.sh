#!/usr/bin/env bash

# Please remember to activate conda environment Covid_Consensus!




################CREDITS#################

#Original pipeline written by Niema Moshiri at UCSD

#This pipeline is a running modification by Sean Sierra-Patev at RIDOH

#Most recent modifications 1-2024

#DISCLAIMER: Tools have been replaced and re-replaced, this pipeline may be a "ship of Thesius". There is no guarantee of continuity.

#This pipeline is named for Conanicut Island in Narragansett Bay. The island's namesake, Canonicus (1565 – June 4, 1647), was Chief Sachem of the Narragansett people at the time of Rhode Island's founding.

########################################


#check directory contents

#FASTA=../ref/*.fasta
#FA=../ref/*.fa
#FNA=../ref/*.fna

if [[ $(ls ../ref/*.{fna,fa,fasta} 2>/dev/null | wc -l) != 1  ]]; then
    echo "Please supply a single FASTA reference in /ref/" ; exit 1
elif [[ $(ls ../ref/*.fna 2>/dev/null | wc -l) == '0' && $(ls ../ref/*.fa 2>/dev/null | wc -l) == '0' && $(ls ../ref/*.fasta 2>/dev/null | wc -l) == '1' ]]; then
  REF_FAS=$(ls ../ref/*.fasta)
elif [[ $(ls ../ref/*.fna 2>/dev/null | wc -l) == '0' && $(ls ../ref/*.fa 2>/dev/null | wc -l) == '1' && $(ls ../ref/*.fasta 2>/dev/null | wc -l) == '0' ]]; then
  REF_FAS=$(ls ../ref/*.fa)
elif [[ $(ls ../ref/*.fna 2>/dev/null | wc -l) == '1' && $(ls ../ref/*.fa 2>/dev/null | wc -l) == '0' && $(ls ../ref/*.fasta 2>/dev/null | wc -l) == '0' ]]; then
  REF_FAS=$(ls ../ref/*.fna)
fi

echo "Provided Reference Genome is" $REF_FAS

#GFF=../ref/*.gff
#GFF3=../ref/*.gff3

if [[ $(ls ../ref/*.{gff,gff3} 2>/dev/null | wc -l) != 1  ]]; then
    echo "Please supply a single gff annotation in /ref/" ; exit 1
elif [[ $(ls ../ref/*.gff 2>/dev/null | wc -l) == '0' && $(ls ../ref/*.gff3 2>/dev/null | wc -l) == '1' ]]; then
  REF_GFF=$(ls ../ref/*.gff3)
elif [[ $(ls ../ref/*.gff 2>/dev/null | wc -l) == '1' && $(ls ../ref/*.gff3 2>/dev/null | wc -l) == '0' ]]; then
  REF_GFF=$(ls ../ref/*.gff)
fi

echo "Provided Reference Assembly is" $REF_GFF

if [[ $(ls ../ref/*.bed 2>/dev/null | wc -l) != 1  ]]; then
    echo "Please supply a single primer BED in /ref/" ; exit 1
fi

PRIMER_BED=$(ls ../ref/*.bed)

echo "Provided Reference Primer Sequences are" $PRIMER_BED

# default constants
THREADS=1
#REF_FAS='../ref/*.fna'
#REF_GFF='../ref/*.gff'
#PRIMER_BED='../ref/*.bed'

# check usage
if [[ "$#" -ne 1 && "$#" -ne 2 ]] ; then
    echo "USAGE: $0 sample_prefix [num_threads]" ; exit 1
elif [[ "$#" -eq 2 ]] ; then
    THREADS=$2
fi

###### New Features #######

# Step 0: fastQC
mkdir ${1}_fastQC
{ time ( fastqc -o ${1}_fastQC -t 2 ${1}_L001_R1_001.fastq.gz ${1}_L001_R2_001.fastq.gz ) ; } 2> ${1}.log.0.fastQC.log

###########################

# Step 1a: Map Reads
{ time ( bbmap.sh -Xmx20g t=$THREADS maxindel=500 ref=$REF_FAS in1=${1}_L001_R1_001.fastq.gz in2=${1}_L001_R2_001.fastq.gz out=${1}.sam ) ; } 2> ${1}.log.1.map.log

# Step 2: Sort Reads
{ time ( samtools sort --threads $THREADS -o ${1}.sorted.bam ${1}.sam) ; } 2> ${1}.log.2.sorted.log

# Step 3: Trim Sorted BAM
{ time ( ivar trim -x 5 -e -i ${1}.sorted.bam -b "$PRIMER_BED" -p ${1}.trimmed ) ; } > ${1}.log.3.trim.log 2>&1

# Step 4: Sort Trimmed BAM
{ time ( samtools sort --threads $THREADS -o ${1}.trimmed.sorted.bam ${1}.trimmed.bam && rm ${1}.trimmed.bam ) ; } 2> ${1}.log.4.sorttrimmed.log

# Step 5: Generate Pile-Up
{ time ( samtools mpileup -A -B -aa -d 0 -Q 0 -f $REF_FAS ${1}.trimmed.sorted.bam ) ; } > ${1}.trimmed.sorted.pileup.txt 2> ${1}.log.5.pileup.log

# Step 6: Call Variants
{ time ( cat ${1}.trimmed.sorted.pileup.txt | ivar variants -r "$REF_FAS" -g "$REF_GFF" -p ${1}.trimmed.sorted.pileup.variants.tsv -m 10 ) ; } 2> ${1}.log.6.variants.log

# Step 7: Call Consensus
{ time ( cat ${1}.trimmed.sorted.pileup.txt | ivar consensus -p ${1}.trimmed.sorted.pileup.consensus -q 20 -m 70 -n N -t 0.5 ) ; } > ${1}.log.7.consensus.log 2>&1

# Step 8: Call Depth
{ time ( samtools depth -d 0 -Q 0 -q 0 -aa ${1}.trimmed.sorted.bam ) ; } > ${1}.sorted.depth.txt 2> ${1}.log.8.depth.log

# Step 9: Qualimap
{ time ( qualimap bamqc -bam ${1}.sorted.bam -nt $THREADS --java-mem-size=4G -outdir ${1}.sorted.stats && tar c ${1}.sorted.stats | pigz -9 -p $THREADS > ${1}.sorted.stats.tar.gz && rm -rf ${1}.sorted.stats ) ; } > ${1}.log.9.qualimap.sorted.log 2>&1

#Step 10: Add Sequence to multifasta
cat ${1}.trimmed.sorted.pileup.consensus.fa >> multi.fasta

# Step 11: Compress
tar -cf ${1}.tar.xz -I pixz $1* --remove-files

