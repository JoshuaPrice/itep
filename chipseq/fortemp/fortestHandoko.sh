#!/bin/bash

## Establish Handoko Control
DIR_CONTROL=/data2/josh/chipseq/controlHandoko
# URL_CONTROL=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/SRR172859/SRR172859.fastq.gz
# mkdir $DIR_CONTROL
# wget -O control_raw.fastq.gz $URL_CONTROL
# mv control_raw.fastq.gz $DIR_CONTROL
# gunzip $DIR_CONTROL/control_raw.fastq.gz
# bowtie /data2/josh/genome/mm8_index -q $DIR_CONTROL/control_raw.fastq  -v 2 -m 1 -3 1 -S 2> $DIR_CONTROL/control_mapped.out > $DIR_CONTROL/control_mapped.sam
# sudo mv macs14_MACS_bedGraph/control/macs14_control_afterfiting_all.bdg.gz $DIR_CONTROL
# sudo mv macs14* $DIR
# sudo rm -r macs14_MACS_bedGraph

# Define arrays for links and final directory names
declare -a linksArray=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/SRR172853/SRR172853.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/SRR172855/SRR172855.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/SRR172850/SRR172850.fastq.gz")

declare -a dirsArray=("/data2/josh/chipseq/CTCF_2"
					"/data2/josh/chipseq/POLR2A"
					"/data2/josh/chipseq/p300")

# get length of the arrays
arraylength=${#linksArray[@]}

# use for loop to read all values and indexes
for (( i=0; i<${arraylength}; i++ ));
do
  echo $i " / " ${arraylength} " : " ${dirsArray[$i]}
    # Make new directory for the files relevant to this TF/HisMod
	DIR=${dirsArray[$i]}
	mkdir $DIR

	# URLs where file to download is located
	URL_EXP=${linksArray[$i]}

	# Download files into /data2
	wget -O exp_raw.fastq.gz $URL_EXP
	mv exp_raw.fastq.gz $DIR

	# Gunzip files in /data2
	gunzip $DIR/exp_raw.fastq.gz

	# Map experiment
	bowtie /data2/josh/genome/mm8_index -q $DIR/exp_raw.fastq  -v 2 -m 1 -3 1 -S 2> $DIR/exp_mapped.out > $DIR/exp_mapped.sam

	# Call peaks
	macs14 -t $DIR/exp_mapped.sam -c $DIR_CONTROL/control_mapped.sam --format SAM  --gsize 2700000000 --name "macs14"  --bw 250 --keep-dup 1 --bdg --single-profile --diag &> $DIR/MACS.out

	sudo mv macs14_MACS_bedGraph/control/macs14_control_afterfitting_all.bdg.gz $DIR
	sudo mv macs14_MACS_bedGraph/treat/macs14_treat_afterfitting_all.bdg.gz $DIR
	sudo mv macs14* $DIR
	sudo rm -r macs14_MACS_bedGraph

	# wc -l $DIR/macs14_peaks.bed
done