# ITEP: Pipeline for Inter-TAD and Intra-TAD Enhancer-Pormoter Analysis
This repo contains the code for my (Joshua Price) honors thesis in Molecular and Cell Biology at UC Berkeley, completed over the course of 2018. This bioinformatics project was conducted in the Tjian-Darzacq group at UC Berkeley under the generous guidance of Prof. Xavier Darzacq and graduate student Maxime Woringer. I wrote all of the code for the project (excluding external libraries) and completed the project from ideation of the research topic to submission of the [thesis paper](thesis.pdf).

## Research Poster
This is the poster I presented at UC Berkeley on November 30, 2018.
![Alt text](poster.jpg?raw=true "ITEP Poster")

## Workflow outline
The data processing procedure is shown graphically in the poster above and described in depth in my [thesis paper](thesis.pdf). Most scripts were written in Python 2.7 with some Linux bash scripts scattered  throughout. For those interested in reading or using the code behind this project, I will summarize the purpose and approach of core files in this repo below. That being said, this project was not designed to be a ready-to-use library for other research groups and many files remain as IPython development notebooks. If your group is interested in replicating or expanding on this work, please reach out to me directly at jprice@berkeley.edu so I can asssit you in converting files as needed for your purposes.

### Step 1: Getting enrichment peaks from raw ChIP-seq reads
See the bash script chipseqrun.sh. This file begins with designated web links that contain the .fastq data files of interest, 
