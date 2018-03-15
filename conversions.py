import pandas as pd
import numpy as np

def convertGeneName(names, currentFormat="NCBI", returnFormat="Ensembl"):
	# Load conversions file
	tfConversionsLoc = '/home/josh/mlep/biopy/data/tfConversions.csv'
	dfConversions = pd.read_csv(tfConversionsLoc)

	# Counters for calculating successful conversion later
	numMatch, totNum = 0, 0 
	if currentFormat.upper() == "NCBI" and returnFormat.upper() == "ENSEMBL":
		convertedNames = pd.DataFrame('', index=np.arange(len(names)), columns=["EnsemblGeneID"])
		for idx, row in names.iterrows():
			ncbiNum = row[0]
			if ncbiNum in dfConversions['NCBI gene ID'].values:
				ensemblGeneName = dfConversions.loc[dfConversions['NCBI gene ID'] == ncbiNum, "Gene stable ID"].values[0]
				convertedNames.loc[idx, "EnsemblGeneID"] = ensemblGeneName
				numMatch += 1
				totNum += 1
			else:
				totNum += 1
		conversionRatio = float(numMatch)/totNum
		print "Fraction of TF from TcoF-DB in Ensembl database: ", conversionRatio
		return convertedNames, conversionRatio
		
	# If not NCBI -> Ensembl ID, return "Error"
	print "Requested conversion not yet implemented!"
	return "Error", -1

# Simple function for testing
def simpleObject():
	ncbiGenes = [14739, 66714, 71911] # NCBI genes
	ncbiDf = pd.DataFrame(data={'NCBI':ncbiGenes})
	return ncbiDf

if __name__ == "__main__":
    a = simpleObject()
    b, r = convertGeneName(a)
    print(b)