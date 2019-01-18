#!/bin/bash

## Establish Wei Control
DIR_CONTROL=/data2/josh/chipseq/controlWei
URL_CONTROL=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001996/SRR001996.fastq.gz
mkdir $DIR_CONTROL
wget -O control_raw.fastq.gz $URL_CONTROL
mv control_raw.fastq.gz $DIR_CONTROL
gunzip $DIR_CONTROL/control_raw.fastq.gz
bowtie /data2/josh/genome/mm8_index -q $DIR_CONTROL/control_raw.fastq  -v 2 -m 1 -3 1 -S 2> $DIR_CONTROL/control_mapped.out > $DIR_CONTROL/control_mapped.sam
sudo mv macs14_MACS_bedGraph/control/macs14_control_afterfiting_all.bdg.gz $DIR_CONTROL
sudo mv macs14* $DIR
sudo rm -r macs14_MACS_bedGraph

# Define arrays for links and final directory names
declare -a linksArray=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002023/SRR002023.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002000/SRR002000.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002012/SRR002012.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002004/SRR002004.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002039/SRR002039.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001992/SRR001992.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001985/SRR001985.fastq.gz"
					"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002027/SRR002027.fastq.gz")

declare -a dirsArray=("/data2/josh/chipseq/Sox2"
					"/data2/josh/chipseq/Klf4"
					"/data2/josh/chipseq/Oct4"
					"/data2/josh/chipseq/Nanog"
					"/data2/josh/chipseq/c-Myc"
					"/data2/josh/chipseq/Esrrb"
					"/data2/josh/chipseq/CTCF"
					"/data2/josh/chipseq/Suz12")

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
	macs14 -t $DIR/exp_mapped.sam -c $DIR_CONTROL/control_mapped.sam --format SAM  --gsize 2700000000 --name "macs14"  --bw 500 --keep-dup 1 --bdg --single-profile --diag &> $DIR/MACS.out

	sudo mv macs14_MACS_bedGraph/control/macs14_control_afterfitting_all.bdg.gz $DIR
	sudo mv macs14_MACS_bedGraph/treat/macs14_treat_afterfitting_all.bdg.gz $DIR
	sudo mv macs14* $DIR
	sudo rm -r macs14_MACS_bedGraph

	# wc -l $DIR/macs14_peaks.bed
done