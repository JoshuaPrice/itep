# ITEP: Pipeline for Inter-TAD and Intra-TAD Enhancer-Pormoter Analysis
This repo contains the code for my (Joshua Price's) honors thesis in Molecular and Cell Biology at UC Berkeley, completed over the course of 2018. This bioinformatics project was conducted in the Tjian-Darzacq group at UC Berkeley under the generous guidance of Prof. Xavier Darzacq and graduate student Maxime Woringer. I wrote all of the code for the project (excluding external libraries) and completed the project from ideation of the research topic to submission of the [thesis paper](thesis.pdf).

## Research Poster
This is the poster I presented at UC Berkeley on November 30, 2018.
![Alt text](poster.jpg?raw=true "ITEP Poster")

## Workflow outline
The data processing procedure is shown graphically in the poster above and described in depth in my [thesis paper](thesis.pdf). Most scripts were written in Python 2.7 with some Linux bash scripts scattered  throughout. For those interested in reading or using the code behind this project, I will summarize the purpose and approach of core files in this repo below. That being said, this project was not designed to be a ready-to-use library for other research groups and many files remain as IPython development notebooks. If your group is interested in replicating or expanding on this work, please reach out to me directly at jprice@berkeley.edu so I can asssit you in converting files as needed for your purposes.

### Step 1: Getting enrichment peaks from raw ChIP-seq reads
See the bash script [chipseqrun.sh](chipseq/chipseqrun.sh). This file begins with designated web links that contain the .fastq data files of interest, maps the experiment reads to the control reference using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and calls peaks using [MACS 14](http://liulab.dfci.harvard.edu/MACS/00README.html). The resulting file .bed file contains peak information, with each row representing a single peak and containing the start and end positions of the enrichment peak as well as the strength of the peak (as determined by MACS). 

Separate files for ChIP-seq experiments that shared a single control are included in the [forchip](chipseq/forchip) directory. Each of these follows the same procedure as in the single experiment approach, but handles multiple experiments by loops through each of the experiments individually.

### Step 2: Establishing an enhancer-promoter interaction (EPI) set for the analysis


### Step 3: Calculating enrichment frequencies of EPIs in the set for transcription factors and histone modifications of interest

--
This summarizes a year's worth of work, so there is naturally much more contained in this repo as well as approaches that were less successful. If you would like to learn more about the project, please refer to my [thesis paper](thesis.pdf) and/or contact me at jprice@berkeley.edu.
