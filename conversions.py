# File name: conversions.py
# Author: Joshua Price
# Date created: 03/08/2018
# Date last modified: 03/15/2018
# Python Version: 2.7

# Purpose: convertGeneName method converts between NCBI, Ensembl, and gene name formats
	# convertGeneName: takes df with gene name, NCBI, or Ensembl format, returns df with converted names
	# ncbiObject: returns df containing some NCBI gene IDs to convert
	# namesObject: returns df containing some gene names to convert

import pandas as pd
import numpy as np

def convertGeneName(names, currentFormat="NCBI", returnFormat="Ensembl"):
	# Purpose: Converts between NCBI ID, Ensembl ID, and gene name formats
	# Inputs: 
		# names: pd df with one column containing NCBI IDs, Ensembl IDs, or gene acronyms
		# currentFormat: indicates format in 'names' - should be "NCBI", "Ensembl", or "Name"
		# returnFormat: indicates format to return - - should be "NCBI", "Ensembl", or "Name"
	# Outputs:
		# 1. convertedNames: pd df with one coulmn containing NCBI IDs, Ensembl IDs, or gene acronyms
		# 2. conversionRatio - Fraction of IDs successfully converted

	# Ensure correct input formatting
	currentFormat = currentFormat.upper()
	returnFormat = returnFormat.upper()
	assert (currentFormat == "NCBI" or currentFormat == "NAME"),"Invalid current format: need either NCBI ID or Gene Name"
	assert (returnFormat == "ENSEMBL"),"Invalid desired format: only Ensembl supported"

	# Import conversions file, originally downloaded from Ensembl
	tfConversionsLoc = '/data2/josh/tfConversions.csv'
	dfConversions = pd.read_csv(tfConversionsLoc)

	# If NCBI format provided and ENSEMBL desired, make that conversion
	numMatch, totNum = 0, 0 
	if currentFormat == "NCBI" and returnFormat == "ENSEMBL":

		# Create new df containing converted names and populate it; to be returned at end
		convertedNames = pd.DataFrame('', index=np.arange(len(names)), columns=["EnsemblGeneID"])
		for idx, row in names.iterrows():
			ncbiNum = row[0] # should only be one row in input
			if ncbiNum in dfConversions['NCBI gene ID'].values:
				ensemblGeneName = dfConversions.loc[dfConversions['NCBI gene ID'] == ncbiNum, "Gene stable ID"].values[0]
				convertedNames.loc[idx, "EnsemblGeneID"] = ensemblGeneName
				numMatch += 1
				totNum += 1
			else:
				totNum += 1
		conversionRatio = float(numMatch)/totNum

		# Report conversion success rate (should be >0.9)
		print "Fraction successfully converted from NCBI to Ensembl: ", conversionRatio
		return convertedNames, conversionRatio
		
	# If name provided and ENSEMBL desired, make that conversion
	elif currentFormat == "NAME" and returnFormat == "ENSEMBL":

		# Create new df containing converted names and populate it; to be returned at end
		convertedNames = pd.DataFrame('', index=np.arange(len(names)), columns=["EnsemblGeneID"])
		for idx, row in names.iterrows():
			geneName = row[0] # should only be one row in input
			if geneName in dfConversions['Gene name'].str.upper().values:
				ensemblGeneName = dfConversions.loc[dfConversions['Gene name'].str.upper() == geneName, "Gene stable ID"].values[0]
				convertedNames.loc[idx, "EnsemblGeneID"] = ensemblGeneName
				numMatch += 1
				totNum += 1
			else:
				totNum += 1
		conversionRatio = float(numMatch)/totNum

		# Report conversion success rate (should be >0.9)
		print "Fraction successfully converted from Name to Ensembl: ", conversionRatio
		return convertedNames, conversionRatio

	# If conversion not supported, throw error
	print "Requested conversion not yet implemented!"
	return "Error", -1

# Simple functions for testing
def ncbiObject():
	ncbiGenes = [14739, 66714, 71911] # NCBI genes
	ncbiDf = pd.DataFrame(data={'NCBI':ncbiGenes})
	return ncbiDf

def namesObject():
	geneNames = ['SOX2', 'KLF4', 'POU5F1']
	geneNameDf = pd.DataFrame(data={'Name':geneNames})
	return geneNameDf


if __name__ == "__main__":
    a = namesObject()
    b, r = convertGeneName(a, "Name", "Ensembl")
    print(b)