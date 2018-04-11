# File name: genLabeledLoops.py
# Author: Joshua Price
# Date created: 04/03/2018
# Date last modified: 04/11/2018
# Python Version: 2.7

# Purpose: addTadCols and addEPCols annotates known contacts with TAD labels and enhancer/promoter labels
	# addTadCols: takes BedTool with TAD information (tads) and df with loops information, adds two TAD cols to loops df
	# addEPCols: takes BedTool with regulatory information (regs) and, df with loops information, adds two reg cols to loops df
	# mergeLabeledLoops: merges two dfs and saves single df with loops, TAD, and reg info into .csv file
	# locFilter: method for pbt filtering, returns true if feature encapsulates given bp info
	# approxLocFilter: method for pbt filtering, returns true if feature overlaps with given bp info
	# containsPromoter: method for pbt filtering, returns true if feature labeled promoter or promoter proximal region
	# containsEnhancer: method for pbt filtering, returns true if feature labeled enhancer
	# extractLoopInfo: helper function that extracts useful contact info from df row
	# bedCopyAndLen: helper function to get BedTool length without erasing it

import pybedtools as pbt
import pandas as pd

def addTadCols(tads, loops):
	# Purpose: takes BedTool with TAD information (tads) and df with loops information, adds two TAD cols to loops df
	# Inputs: 
		# tads: BedTool with TAD information in bed format (chr#, startBp, endBp, name, score, ...)
		# loops: pd df (no header) with chrA#, startBpA, endBpA, chrB#, startBpB, endBpB, score
	# Outputs:
		# loops: pd df with two new columns (tadA, tadB) populated with TAD IDs corresponding to loops

	# Add new columns to loops df
	loops["tadA"] = ""
	loops["tadB"] = ""

	# Initialize counters used to produce early report
	numLoops, numIntra, numProblem = 0, 0, 0 

	for idx, row in loops.iterrows():
		chrA, startBpA, endBpA, chrB, startBpB, endBpB = extractLoopInfo(row)

		# tadAInfo contains TAD(s) that encapsulate the contact point
		tadAInfo = tads.filter(locFilter, chrA, startBpA, endBpA)
		tadAInfo, tadAInfoLen = bedCopyAndLen(tadAInfo)

		# if contact point encapsulated by TAD, add that TAD label to df
		if tadAInfoLen > 0: 
			loops.loc[idx, "tadA"] = tadAInfo[0].name

		# tadBInfo contains TAD(s) that encapsulate the contact point
		tadBInfo = tads.filter(locFilter, chrB, startBpB, endBpB)
		tadBInfo, tadBInfoLen = bedCopyAndLen(tadBInfo)

		# if contact point encapsulated by TAD, add that TAD label to df
		if tadBInfoLen > 0: 
			loops.loc[idx, "tadB"] = tadBInfo[0].name

		# if both contact points in same TAD (and in a TAD at all), note as intra-TAD
		if loops.loc[idx, "tadA"] == loops.loc[idx, "tadB"] and tadAInfoLen > 0:
			numIntra += 1
		
		# if either contact point is not within a TAD, note as TADless ('problem') contact
		if tadAInfoLen == 0 or tadBInfoLen == 0:
			numProblem += 1

		numLoops += 1
		print "Percent intra-TAD so far: ", float(numIntra)/(numLoops-numProblem)
		print "Rows completed: ", numLoops

	# Provide initial report about fraction of intra-TAD contacts and fraction of TADless contacts
	print "Percent intra-TAD: ", float(numIntra)/(numLoops-numProblem)
	print "Percent overlapping with TAD boundary (not counted in inter/intra): ", float(numProblem) / numLoops
	return loops

def addEPCols(regs, loops):
	# Purpose: takes BedTool with regulatory information (regs) and df with loops information, adds two reg cols to loops df
	# Inputs: 
		# regs: BedTool with regulatory information in bed format (chr#, startBp, endBp, name)
		# loops: pd df (no header) with chrA#, startBpA, endBpA, chrB#, startBpB, endBpB, score, (tadA), (tadB)
	# Outputs:
		# loops: pd df with two new columns (regA, regB) populated with regulatory feature IDs corresponding to loops

	# Add new columns to loops df	
	loops["regA"] = ""
	loops["regB"] = ""

	# Initialize counters used to produce early report
	numLoops, numEE, numPP, numEP, numProblem = 0, 0, 0, 0, 0

	for idx, row in loops.iterrows():
		chrA, startBpA, endBpA, chrB, startBpB, endBpB = extractLoopInfo(row)

		# regAInfo contains regulatory labels of regions that overlap with the contact point
		regAInfo = regs.filter(approxLocFilter, chrA, startBpA, endBpA)
		regAInfo, regAInfoLen = bedCopyAndLen(regAInfo)

		# labels for whether a promoter or enhancer labeled region overlaps with the contact point
		promoterPresentA = len(regAInfo.filter(containsPromoter)) > 0
		enhancerPresentA = len(regAInfo.filter(containsEnhancer)) > 0

		# add regulatory label to df as appropriate (EP = both enhancer AND promoter, N = neither)
		if promoterPresentA and enhancerPresentA:
			loops.loc[idx, "regA"] = 'EP'
		elif promoterPresentA:
			loops.loc[idx, "regA"] = 'P'
		elif enhancerPresentA:
			loops.loc[idx, "regA"] = 'E'
		else:
			loops.loc[idx, "regA"] = 'N'

		# regBInfo contains regulatory labels of regions that overlap with the contact point
		regBInfo = regs.filter(approxLocFilter, chrB, startBpB, endBpB)
		regBInfo, regBInfoLen = bedCopyAndLen(regBInfo)

		# labels for whether a promoter or enhancer labeled region overlaps with the contact point
		promoterPresentB = len(regBInfo.filter(containsPromoter)) > 0
		enhancerPresentB = len(regBInfo.filter(containsEnhancer)) > 0

		# add regulatory label to df as appropriate (EP = both enhancer AND promoter, N = neither)
		if promoterPresentB and enhancerPresentB:
			loops.loc[idx, "regB"] = 'EP'
		elif promoterPresentB:
			loops.loc[idx, "regB"] = 'P'
		elif enhancerPresentB:
			loops.loc[idx, "regB"] = 'E'
		else:
			loops.loc[idx, "regB"] = 'N'

		# Count contacts for initial report
		if (promoterPresentA and enhancerPresentB) or (enhancerPresentA and promoterPresentB):
			numEP += 1
		if enhancerPresentA and enhancerPresentB:
			numEE += 1
		if promoterPresentA and promoterPresentB:
			numPP += 1
		numLoops += 1

		# Provide initial report about regulatory distribution of loops
		print "Percent EE: ", float(numEE)/(numLoops)
		print "Percent PP: ", float(numPP)/(numLoops)
		print "Percent EP: ", float(numEP)/(numLoops)
		print "Rows completed: ", numLoops

	# returned original loops df with cols (regA and regB) added
	return loops

def mergeLabeledLoops(tadLabelsDf, regLabelsDf):
	# Purpose: takes BedTool with regulatory information (regs) and df with loops information, adds two reg cols to loops df
	# Inputs: 
		# regs: BedTool with regulatory information in bed format (chr#, startBp, endBp, name)
		# loops: pd df (no header) with chrA#, startBpA, endBpA, chrB#, startBpB, endBpB, score, (tadA), (tadB)
	# Outputs:
		# loops: pd df with two new columns (regA, regB) populated with regulatory feature IDs corresponding to loops

	wholeLabeledDf = pd.concat([tadLabelsDf, regLabelsDf['regA'], \
		regLabelsDf['regB']], axis=1)
	# wholeLabeledDf.to_csv('/data2/josh/expCH12/labeled_loops.csv', index=True, sep='\t', header=True)
	return wholeLabeledDf

def locFilter(feature, chrom, startBp, endBp):
	# pybedtools filter function: Returns True if feature overlaps with contact
	return str(feature.chrom) == chrom and feature.start < startBp and feature.end > endBp

def approxLocFilter(feature, chrom, startBp, endBp):
	# pybedtools filter function: Returns True if feature overlaps with contact
	chrom = 'chr' + str(chrom)
	return str(feature.chrom) == chrom and ((feature.start < startBp and feature.end > startBp) \
		or (feature.start < endBp and feature.end > endBp)
		or (feature.start > startBp and feature.end < endBp))

def containsPromoter(feature):
	# pybedtools filter function: Returns True if feature is promoter
	return feature.name == 'feature_type=Promoter'

def containsEnhancer(feature):
	# pybedtools filter function: Returns True if feature is enhancer
	return feature.name == 'feature_type=Enhancer'

def extractLoopInfo(row):
	# helper function: extracts bed-formatted information from pandas df row
	chrA = row["chrA"]
	chrA = chrA[3:] # remove 'chr'
	startBpA = row["startbpA"]
	endBpA = row["endbpA"]
	chrB = row["chrB"]
	chrB = chrB[3:] # remove 'chr'
	startBpB = row["startbpB"]
	endBpB = row["endbpB"]
	return chrA, startBpA, endBpA, chrB, startBpB, endBpB

def bedCopyAndLen(bTool):
	# helper function: copies bedTool to report length (lame but necessary with pybedtools)
	bToolCopy = bTool.saveas('btool.bed')
	bToolLen = len(bToolCopy)
	bTool = pbt.BedTool('btool.bed')
	return bTool, bToolLen


if __name__ == "__main__":
	# Use CH12 TADs BED file for now
	tadFile = '/data2/josh/expCH12/CH12_lieberman_intra_5kb_domains.bed'
	tads = pbt.BedTool(tadFile)

	loopFile = '/data2/josh/expCH12/lieberman_loops_CH12.txt'
	loops = pd.read_csv(loopFile, sep="\t", header=None)
	loops.columns = ["chrA", "startbpA", "endbpA", "chrB", "startbpB", "endbpB", "score"]

	tadLabeledLoops = addTadCols(tads, loops)
	tadLabeledLoops.to_csv('/data2/josh/expCH12/tad_labeled_loops.csv', index=True, sep='\t', header=True)

	regFile = '/data2/josh/expCH12/mm9_regulatory_converted.bed'
	regs = pbt.BedTool(regFile)
	addEPCols(regs, tadLabeledLoops)
	loops.to_csv('/data2/josh/expCH12/labeled_loops.csv', index=True, sep='\t', header=True)
