from pandas import DataFrame, read_csv
import pandas as pd
import conversions as cnvr

## Part 1: Read .xlsx file

# Path to Excel file with TF names and IDs
# dfTfList = full TF list from TcoF-DB v2: http://tools.sschmeier.com/tcof/home/
tfListLoc = '/home/josh/mlep/biopy/data/allTF.xlsx'
dfTfList = pd.read_excel(tfListLoc)
dfTfList = dfTfList.drop([u'Name', u'Species', u'Type'], axis=1)
del dfTfList[u'# PPI\xa0'] # had difficulty using drop method on this name

# Add ensembl name to dfTfList
ncbiDf = DataFrame(dfTfList["GeneID"])
dfTfList["EnsemblGeneId"], cr = cnvr.convertGeneName(ncbiDf, "NCBI", "Ensembl")

# Add column with RTCC Value
rnaSeqLoc = '/home/josh/mlep/biopy/data/mES_RNAseq.xlsx'
dfRnaSeq = pd.read_excel(rnaSeqLoc, 0)
dfTfList["RTCC"] = ""

numMatchRna = 0
totNumRna = 0
for idx, row in dfTfList.iterrows():
	ensemblNum = row["EnsemblGeneId"]
	if ensemblNum in dfRnaSeq["EnsemblGeneId"].values:
		numMatchRna += 1
		totNumRna += 1
		rtccVal = dfRnaSeq.loc[dfRnaSeq["EnsemblGeneId"] == ensemblNum, "RTCC07_value_1"].values[0]
		dfTfList.loc[idx, "RTCC"] = rtccVal
	else:
		# If row has no RNA-seq peak, give it a 0
		dfTfList.loc[idx, "RTCC"] = 0 # didn't show up in RNA-seq
		totNumRna += 1

print "Fraction of TF from TcoF-DB+Ensembl in RNA-seq: ", float(numMatchRna)/totNumRna

dfTfList.set_index(dfTfList["EnsemblGeneId"], inplace=True)
dfTfList = dfTfList.drop(['EnsemblGeneId'], axis=1)

sortedTf = dfTfList.sort_values(by=['RTCC'], ascending=False)



# Save to csv file
sortedTf.to_csv('data/sortedTF.csv', index=True, header=True)
