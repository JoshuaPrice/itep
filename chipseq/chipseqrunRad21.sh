#!/bin/bash
########################
#
# Chip-Seq Script
# Written by: Joshua Price, jprice@berkeley.edu
# Initiated on: 10/16/18
# Last edited: 11/27/18
#
########################

# Make new directory for the files relevant to this TF/HisMod
DIR=/data2/josh/chipseq/Rad21
mkdir $DIR

# URLs where file to download is located
URL_EXP=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR331/009/SRR3313239/SRR3313239.fastq.gz

# URL for control
URL_CONTROL=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR331/000/SRR3313240/SRR3313240.fastq.gz

# Download files into /data2
wget -O exp_raw.fastq.gz $URL_EXP
wget -O control_raw.fastq.gz $URL_CONTROL

mv exp_raw.fastq.gz $DIR
mv control_raw.fastq.gz $DIR

# Gunzip files in /data2
gunzip $DIR/exp_raw.fastq.gz
gunzip $DIR/control_raw.fastq.gz

# Prepare index file - already done
# bowtie-build /data2/josh/genome/Mus_musculus.NCBIM37.67.dna.toplevel.fa /data2/josh/genome/mm9_index

# Map experiment and control
bowtie /data2/josh/genome/mm9_index -q $DIR/exp_raw.fastq  -v 2 -m 1 -3 1 -S 2> $DIR/exp_mapped.out > $DIR/exp_mapped.sam
bowtie /data2/josh/genome/mm9_index -q $DIR/control_raw.fastq  -v 2 -m 1 -3 1 -S 2> $DIR/control_mapped.out > $DIR/control_mapped.sam

# Call peaks
macs14 -t $DIR/exp_mapped.sam -c $DIR/control_mapped.sam --format SAM  --gsize 2700000000 --name "macs14"  --bw 250 --keep-dup 1 --bdg --single-profile --diag &> $DIR/MACS.out

sudo mv macs14_MACS_bedGraph/control/macs14_control_afterfiting_all.bdg.gz $DIR
sudo mv macs14_MACS_bedGraph/treat/macs14_treat_afterfiting_all.bdg.gz $DIR
sudo mv macs14* $DIR
sudo rm -r macs14_MACS_bedGraph

# wc -l $DIR/macs14_peaks.bed
