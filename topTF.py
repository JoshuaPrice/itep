from pandas import DataFrame
import pandas as pd

# Find expression values of top transcription factors
topTfNames = ['Sox2', 'Klf4', 'Pou5f1']
topTfList = ['ENSMUSG00000074637', 'ENSMUSG00000003032', 'ENSMUSG00000024406']
topDf = DataFrame(data={'EnsemblGene':topTfList, 'Names':topTfNames})
topDf.set_index(topDf["EnsemblGene"], inplace=True)
topDf = topDf.drop([u'EnsemblGene'], axis=1)

sortedTfLoc = '/home/josh/mlep/biopy/data/sortedTF.csv'
sortedTf = pd.read_csv(sortedTfLoc)
sortedTf.set_index(sortedTf["EnsemblGene"], inplace=True)
sortedTf = sortedTf.drop(['EnsemblGene'], axis=1)

topDf["RTCC"] = ""

numMatch = 0
totNum = 0
for ensemblNum, row in topDf.iterrows():
	if ensemblNum in sortedTf.index:
		numMatch += 1
		totNum += 1
		peakVal = sortedTf.loc[sortedTf.index == ensemblNum, "RTCC"].values[0]
		topDf.loc[ensemblNum, "RTCC"] = peakVal
	else:
		# If row has no RNA-seq peak, give it a 0
		topDf.loc[ensemblNum, "RTCC"] = 0 # didn't show up in RNA-seq
		totNum += 1

print(topDf)

print "Fraction of top TF in sorted TF set: ", float(numMatch)/totNum
