# File name: generateTfList.py
# Author: Joshua Price
# Date created: 03/08/2018
# Date last modified: 03/15/2018
# Python Version: 2.7

# Purpose: Generates .csv file containing TF list ranked by RTCC RNA-seq score
	# Part 1: Import TF info and RNA-seq .xlsx files to df's and format df's
	# Part 2: Add column with RTCC score to TF df
	# Part 3: Sort TFs according to RTCC score, save to .csv

import pandas as pd
import conversions as cnvr


### Part 1: Import TF info and RNA-seq .xlsx files to df's and format df's ###

# Full TF list from TcoF-DB v2: http://tools.sschmeier.com/tcof/home/
# Import TF list and format df
tfListLoc = '/home/josh/mlep/data/allTF.xlsx'
tfListDf = pd.read_excel(tfListLoc)
tfListDf = tfListDf.drop([u'Name', u'Species', u'Type'], axis=1)
del tfListDf[u'# PPI\xa0'] # had difficulty using drop method on this name

# Append additional proteins of interest to TF list (such as Polr2a)
otherProteinsLoc = '/home/josh/mlep/data/moreProteinsInterest.xlsx'
otherProteinsDf = pd.read_excel(otherProteinsLoc)
tfListDf = tfListDf.append(otherProteinsDf, ignore_index=True)

# Make sure gene name/acronyms all in upper format
tfListDf['Symbol'] = tfListDf['Symbol'].str.upper()

# Import RNA-seq data, obtained from Claudia
rnaSeqLoc = '/home/josh/mlep/data/mES_RNAseq.xlsx'
rnaSeqDf = pd.read_excel(rnaSeqLoc, 0)

# Add Ensembl ID to tfListDf
ncbiDf = pd.DataFrame(tfListDf["GeneID"])
tfListDf["EnsemblGeneId"], cr = cnvr.convertGeneName(ncbiDf, "NCBI", "Ensembl")


### Part 2: Add and populate column with RTCC score to TF df ###

# Add column with RTCC Value to df
tfListDf["RTCC"] = ""

# Populate RTCC score column by matching Ensembl ID's between df's
numMatchRna, totNumRna = 0, 0
for idx, row in tfListDf.iterrows():
	ensemblNum = row["EnsemblGeneId"]
	if ensemblNum in rnaSeqDf["EnsemblGeneId"].values:
		numMatchRna += 1
		totNumRna += 1
		rtccVal = rnaSeqDf.loc[rnaSeqDf["EnsemblGeneId"] == ensemblNum, "RTCC07_value_1"].values[0]
		tfListDf.loc[idx, "RTCC"] = rtccVal
	else:
		# If row has no RNA-seq peak, give it a 0
		tfListDf.loc[idx, "RTCC"] = float('nan') # didn't show up in RNA-seq
		totNumRna += 1

# Report successful conversion rate (should be > 0.95)
print "Fraction of TF from TcoF-DB+Ensembl in RNA-seq: ", float(numMatchRna)/totNumRna


### Part 3: Sort TFs according to RTCC score, save to .csv ###

sortedTfDf = tfListDf.sort_values(by=['RTCC'], ascending=False)
sortedTfDf.reset_index(drop=True, inplace=True)
sortedTfDf.index.name = "Rank0"

# Save to sortedTf.csv
sortedTfDf.to_csv('data/sortedTf.csv', index=True, header=True)
