from pandas import DataFrame, read_csv
import pandas as pd

## Part 1: Read .xlsx file

# Path to Excel file with TF names and IDs
# dfTfList = full TF list from TcoF-DB v2: http://tools.sschmeier.com/tcof/home/
tfListLoc = '/home/josh/mlep/biopy/data/allTF.xlsx'
dfTfList = pd.read_excel(tfListLoc, 0, index_col='GeneID')
dfTfList = dfTfList.drop([u'Name', u'Species', u'Type'], axis=1)
del dfTfList[u'# PPI\xa0'] # had difficulty using drop method on this name

# Ensembl -> NCBI conversions from Ensembl mart view
tfConversionsLoc = '/home/josh/mlep/biopy/data/tfConversions.csv'
dfConversions = pd.read_csv(tfConversionsLoc, index_col='NCBI gene ID')

# Add ensembl name to dfTfList
dfTfList["EnsemblGene"] = ""

# Add ensembl names to rows where matching NCBI IDs
# and determine percent matching
numMatch = 0
totNum = 0

for ncbiNum, row in dfTfList.iterrows():
	if ncbiNum in dfConversions.index:
		numMatch += 1
		totNum += 1
		ensemblGeneName = dfConversions.loc[dfConversions.index == ncbiNum, "Gene stable ID"].values[0]
		dfTfList.loc[ncbiNum, "EnsemblGene"] = ensemblGeneName
	else:
		totNum += 1

print "Fraction of TF from TcoF-DB in Ensembl database: ", float(numMatch)/totNum

# Change index of dfTfList to EnsemblNum
dfTfList = dfTfList.reset_index()
dfTfList.rename(columns={'index':'NCBI'}, inplace=True)
dfTfList.set_index(dfTfList["EnsemblGene"], inplace=True)
dfTfList = dfTfList.drop(['EnsemblGene'], axis=1)

# Add column with RTCC Value
rnaSeqLoc = '/home/josh/mlep/biopy/data/mES_RNAseq.xlsx'
dfRnaSeq = pd.read_excel(rnaSeqLoc, 0)
dfRna = dfRnaSeq.set_index(dfRnaSeq["EnsemblGeneId"])

dfTfList["RTCC"] = ""

numMatchRna = 0
totNumRna = 0
for ensemblNum, row in dfTfList.iterrows():
	if ensemblNum in dfRna.index:
		numMatchRna += 1
		totNumRna += 1
		peakVal = dfRna.loc[dfRna.index == ensemblNum, "RTCC07_value_1"].values[0]
		dfTfList.loc[ensemblNum, "RTCC"] = peakVal
	else:
		# If row has no RNA-seq peak, give it a 0
		dfTfList.loc[ensemblNum, "RTCC"] = 0 # didn't show up in RNA-seq
		totNumRna += 1

print "Fraction of TF from TcoF-DB+Ensembl in RNA-seq: ", float(numMatchRna)/totNumRna
dfTfList.set_index(dfTfList["EnsemblGene"], inplace=True)
dfTfList = dfTfList.drop(['EnsemblGene'], axis=1)

sortedTf = dfTfList.sort_values(by=['RTCC'], ascending=False)

# Save to csv file
sortedTf.to_csv('data/sortedTF.csv', index=True, header=True)
