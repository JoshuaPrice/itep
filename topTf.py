# File name: topTf.py
# Author: Joshua Price
# Date created: 03/08/2018
# Date last modified: 03/15/2018
# Python Version: 2.7

# Purpose: Displays RTCC scores and "ranks" of pre-defined TF's of interest

import pandas as pd
import conversions as cnvr

# Format pre-determined top TF into df
topTfNames = ['CTCF', 'YY1', 'NSD2', 'WDR5', 'SUZ12', 'RAD21', 'POLR2A', 'MXI1', \
'MAZ', 'JUN', 'GATA2', 'FOXM1', 'FOS', 'CREBBP', 'CHD7', 'CHD4', 'CEBPB', 
'H2AZ', 'SOX2', 'KLF4', 'POU5F1', 'MYC', 'ESRRB', 'FOXA2']
topTfDf = pd.DataFrame(data={'Names':topTfNames})

# Create Ensembl Gene ID row; cr = conversion ratio (not used)
topTfDf["EnsemblGeneId"], cr = cnvr.convertGeneName(topTfDf, "Name", "Ensembl")

# Add columns to df to store RTCC Values and Rank
sortedTfLoc = '/home/josh/mlep/data/sortedTf.csv'
dfSortedTf = pd.read_csv(sortedTfLoc)
topTfDf["RTCC"] = ""
topTfDf["Rank"] = ""

# Populate RTCC score column by matching Ensembl ID's between df's
numMatchRna, totNumRna = 0, 0
for idx, row in topTfDf.iterrows():
	ensemblNum = row["EnsemblGeneId"]
	if ensemblNum in dfSortedTf["EnsemblGeneId"].values:
		numMatchRna += 1
		totNumRna += 1
		rtccVal = dfSortedTf.loc[dfSortedTf["EnsemblGeneId"] == ensemblNum, "RTCC"].values[0]
		topTfDf.loc[idx, "RTCC"] = rtccVal
		tfRank0 = dfSortedTf.loc[dfSortedTf["EnsemblGeneId"] == ensemblNum, "Rank0"].values[0]
		topTfDf.loc[idx, "Rank"] = tfRank0 + 1  # +1 because using 0-index value for tfRank0
	else:
		# If row has no RNA-seq peak, give score of NaN
		topTfDf.loc[idx, "RTCC"] = float('nan')
		topTfDf.loc[idx, "Rank"] = "inf"
		totNumRna += 1

# Report conversion success rate (should be >0.9)
print "Fraction of TF from TcoF-DB+Ensembl in RNA-seq: ", float(numMatchRna)/totNumRna

# Display the RTCC scores and "ranks" of TF's of interest
print topTfDf
