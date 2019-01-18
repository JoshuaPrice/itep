import pybedtools as pbt
from pybedtools import featurefuncs
import pandas as pd
import gffutils
from gffutils import pybedtools_integration

def simpleBedExample():
	a = pbt.example_bedtool('a.bed')
	b = pbt.example_bedtool('b.bed')
	print a.intersect(b)

def printBedFile():
	# Use CH12 TADs BEDfile for now
	tadFile = '/data2/josh/expCH12/CH12_lieberman_intra_5kb_domains.bed'
	a = pbt.BedTool(tadFile)
	# for interval in a:
	# 	print interval
	print a[0]

def simpleFilter():
	# Use CH12 TADs BED file for now
	tadFile = '/data2/josh/expCH12/CH12_lieberman_intra_5kb_domains.bed'
	tads = pbt.BedTool(tadFile)
	print tads

	loopFile = '/data2/josh/expCH12/lieberman_loops_CH12.txt'
	loops = pd.read_csv(loopFile, sep="\t", header=None)
	loops.columns = ["chrA", "startbpA", "endbpA", "chrB", "startbpB", "endbpB", "score"]
	# print loops

	# Add new columns to df
	loops["tadA"] = ""
	loops["tadB"] = ""

	numLoops, numIntra, numProblem = 0, 0, 0 

	for idx, row in loops.iterrows():
		chrA = row["chrA"]
		chrA = chrA[3:]
		startBpA = row["startbpA"]
		endBpA = row["endbpA"]
		tadAInfo = tads.filter(locFilter, chrA, startBpA, endBpA)

		tadAInfoCopy = tadAInfo.saveas('tadAtemp.bed')
		tadAInfoLen = len(tadAInfoCopy)
		tadAInfo = pbt.BedTool('tadAtemp.bed')

		if tadAInfoLen > 0: 
			row["tadA"] = tadAInfo[0].name

		chrB = row["chrB"]
		chrB = chrB[3:]
		startBpB = row["startbpB"]
		endBpB = row["endbpB"]
		tadBInfo = tads.filter(locFilter, chrB, startBpB, endBpB)

		tadBInfoCopy = tadBInfo.saveas('tadBtemp.bed')
		tadBInfoLen = len(tadBInfoCopy)
		tadBInfo = pbt.BedTool('tadBtemp.bed')

		if tadBInfoLen > 0: 
			row["tadB"] = tadBInfo[0].name

		if row["tadA"] == row["tadB"] and tadAInfoLen > 0:
			numIntra += 1

		if tadAInfoLen == 0 or tadBInfoLen == 0:
			numProblem += 1

		numLoops += 1
		print numIntra, numProblem, numLoops
		return  # only first row


def addTadCols():
	# Use CH12 TADs BED file for now
	tadFile = '/data2/josh/expCH12/CH12_lieberman_intra_5kb_domains.bed'
	tads = pbt.BedTool(tadFile)

	loopFile = '/data2/josh/expCH12/lieberman_loops_CH12.txt'
	loops = pd.read_csv(loopFile, sep="\t", header=None)
	loops.columns = ["chrA", "startbpA", "endbpA", "chrB", "startbpB", "endbpB", "score"]
	# print loops

	# Add new columns to df
	loops["tadA"] = ""
	loops["tadB"] = ""

	numLoops, numIntra, numProblem = 0, 0, 0 

	for idx, row in loops.iterrows():
		chrA = row["chrA"]
		chrA = chrA[3:]
		startBpA = row["startbpA"]
		endBpA = row["endbpA"]
		tadAInfo = tads.filter(locFilter, chrA, startBpA, endBpA)

		tadAInfoCopy = tadAInfo.saveas('tadAtemp.bed')
		tadAInfoLen = len(tadAInfoCopy)
		tadAInfo = pbt.BedTool('tadAtemp.bed')

		if tadAInfoLen > 0: 
			loops.loc[idx, "tadA"] = tadAInfo[0].name

		chrB = row["chrB"]
		chrB = chrB[3:]
		startBpB = row["startbpB"]
		endBpB = row["endbpB"]
		tadBInfo = tads.filter(locFilter, chrB, startBpB, endBpB)

		tadBInfoCopy = tadBInfo.saveas('tadBtemp.bed')
		tadBInfoLen = len(tadBInfoCopy)
		tadBInfo = pbt.BedTool('tadBtemp.bed')

		if tadBInfoLen > 0: 
			loops.loc[idx, "tadB"] = tadBInfo[0].name

		if row["tadA"] == row["tadB"] and tadAInfoLen > 0:
			numIntra += 1

		if tadAInfoLen == 0 or tadBInfoLen == 0:
			numProblem += 1

		numLoops += 1
		print numIntra, numProblem, numLoops

	print "Percent intra-TAD: ", float(numIntra)/(numLoops-numProblem)
	print "Percent overlapping with TAD boundary (not counted in inter/intra): ", float(numProblem) / numLoops

	# Save to sortedTf.csv
	loops.to_csv('/data2/josh/expCH12/tad_labeled_loops.csv', index=True, sep='\t', header=True)
	return loops


def locFilter(feature, chrom, startBp, endBp):
	"Returns True if feature TAD encapsulates one side of contact"
	return str(feature.chrom) == chrom and feature.start < startBp and feature.end > endBp

def convertGffToBed(gffFile): 
	# Doesn't work as desired because doesn't change ordering of columns
	toConvert = pbt.BedTool(gffFile)
	toConvert.saveas(gffFile[:-4] + ".bed")


def convertGffToBedPandas(gffFile):
	# Gff features column makes this difficult
	toConvert = pd.read_csv(gffFile, sep="\t", header=None)
	toConvert.columns = ["chr", "source", "feature", "start", "end", "score", \
		"strand", "frame", "attribute"]
	toConvert = toConvert.drop([u'source', u'feature', u'score', \
		'strand', 'frame'], axis=1)

	# add 'chr' to beginning of chromosomes
	toConvert['chr'] = 'chr' + toConvert['chr'].astype(str)

	splitGffAttributes = lambda attributesCol: pd.Series([i for i in reversed(attributesCol.split(';'))])
	attributesSeries = toConvert['attribute'].apply(splitGffAttributes)
	attributesSeries.rename(columns={0:'feature_type'},inplace=True)

	toConvert = toConvert.drop([u'attribute'], axis=1)

	return pd.concat([toConvert, attributesSeries['feature_type']], axis=1)


	

def convertGffToBedGffUtils(gffFile):
	fn = gffutils.example_filename(gffFile)
	# db = gffutils.create_db(fn, dbfn=gffFile[:-4] + '.db', force=True, keep_order=True, \
	# 	merge_strategy='merge', sort_attribute_values=True)
	db = gffutils.FeatureDB(gffFile[:-4] + '.db', keep_order=True)
	gffIterator = db.all_features(order_by='start')
	bedVersion = pybedtools_integration.to_bedtool(gffIterator)
	print bedVersion

	# bedVersion.saveas(gffFile[:-4] + '.bed')

def convertGffToBedPbt(gffFile):
	toConvert = pbt.BedTool(gffFile)
	featureToConvert = toConvert[1:4]
	print "featureToConvert type", type(featureToConvert)
	# bedVersion = featurefuncs.gff2bed(featureToConvert, name_field='ID')
	bedVersion = pbt.BedTool()
	convertedFeatures = featurefuncs.gff2bed(toConvert[0], name_field='feature_type')

	for interval in featureToConvert:
		# print interval
		newInterval = featurefuncs.gff2bed(interval, name_field='feature_type')
		print type(newInterval)
		bedVersion.cat(newInterval)
		print bedVersion
	# print bedVersion
	
	# return



	# for region in db.features_of_type('regulatory_region', order_by='start'):
	#	print region['ID'][0]


if __name__ == "__main__":
    # addTadCols()
    # simpleFilter()
    bedReg = convertGffToBedPandas('/data2/josh/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20161111.gff')
    bedReg.to_csv('/data2/josh/expCH12/mm10_regulatory.bed', index=False, sep='\t', header=False)