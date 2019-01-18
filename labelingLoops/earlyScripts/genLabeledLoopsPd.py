# File name: genLabeledLoopsPd.py
# Author: Joshua Price
# Date created: 04/03/2018
# Date last modified: 04/12/2018
# Python Version: 2.7

# Purpose: addTadCols and addEPCols annotates known contacts with TAD labels and enhancer/promoter labels
	# addTadCols: takes pd df with TAD information (tads) and df with loops information, adds two TAD cols to loops df
	# addEPCols: takes pd df with regulatory information (regs) and, df with loops information, adds three reg cols to loops df

import pybedtools as pbt
import pandas as pd
import numpy as np

def addTadCols(tads, loops):
	# Purpose: takes pd df with TAD information (tads) and df with loops information, adds two TAD cols to loops df
	# Inputs: 
		# tads: pd df with TAD information in bed format (chr#, startBp, endBp, name, score, ...)
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
		tadAInfo = tads[(tads['chr'] == chrA) & (tads['startbp'] < startBpA) & (tads['endbp'] > endBpA)]
		tadAInfoLen = tadAInfo.shape[0]
		tadAInfo.reset_index(drop=True, inplace=True)

		# if contact point encapsulated by TAD, add that TAD label to df
		if tadAInfoLen > 0: 
			loops.loc[idx, "tadA"] = tadAInfo.loc[0, "name"]
		else: 
			loops.loc[idx, "tadA"] = np.NaN

		# tadBInfo contains TAD(s) that encapsulate the contact point
		tadBInfo = tads[(tads['chr'] == chrB) & (tads['startbp'] < startBpB) & (tads['endbp'] > endBpB)]
		tadBInfoLen = tadBInfo.shape[0]
		tadBInfo.reset_index(drop=True, inplace=True)

		# if contact point encapsulated by TAD, add that TAD label to df
		if tadBInfoLen > 0: 
			loops.loc[idx, "tadB"] = tadBInfo.loc[0, "name"]
		else: 
			loops.loc[idx, "tadB"] = np.NaN

		# if both contact points in same TAD (and in a TAD at all), note as intra-TAD
		if loops.loc[idx, "tadA"] == loops.loc[idx, "tadB"] and tadAInfoLen > 0:
			numIntra += 1
		
		# if either contact point is not within a TAD, note as TADless ('problem') contact
		if tadAInfoLen == 0 or tadBInfoLen == 0:
			numProblem += 1

		numLoops += 1
		# print "Percent intra-TAD so far: ", float(numIntra)/(numLoops-numProblem)
		# print "Rows completed: ", numLoops

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
	loops["Int_Status"] = ""

	# Initialize counters used to produce early report
	numLoops, numEE, numPP, numEP, numProblem = 0, 0, 0, 0, 0

	for idx, row in loops.iterrows():
		chrA, startBpA, endBpA, chrB, startBpB, endBpB = extractLoopInfo(row)
		chrA = 'chr' + str(chrA)
		chrB = 'chr' + str(chrB)

		# regAInfo contains regulatory labels of regions that overlap with the contact point        
		regAInfo = regs[(((regs['chr'] == chrA) & (regs['startbp'] < startBpA) \
				& (regs['endbp'] > startBpA)) | \
					((regs['chr'] == chrA) & (regs['startbp'] < endBpA) & (regs['endbp'] > endBpA)) | \
					((regs['chr'] == chrA) & (regs['startbp'] > startBpA) & (regs['endbp'] < endBpA)))]

		regAInfoLen = regAInfo.shape[0]
		regAInfo.reset_index(drop=True, inplace=True)

		# labels for whether a promoter or enhancer labeled region overlaps with the contact point
		promoterPresentA = regAInfo[regAInfo['feature_type'] == 'feature_type=Promoter'].shape[0] > 0
		enhancerPresentA = regAInfo[regAInfo['feature_type'] == 'feature_type=Enhancer'].shape[0] > 0
        
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
		regBInfo = regs[(((regs['chr'] == chrB) & (regs['startbp'] < startBpB) \
				& (regs['endbp'] > startBpB)) | \
					((regs['chr'] == chrB) & (regs['startbp'] < endBpB) & (regs['endbp'] > endBpB)) | \
					((regs['chr'] == chrB) & (regs['startbp'] > startBpB) & (regs['endbp'] < endBpB)))]
        
		regBInfoLen = regBInfo.shape[0]
		regBInfo.reset_index(drop=True, inplace=True)

		# labels for whether a promoter or enhancer labeled region overlaps with the contact point
		promoterPresentB = regBInfo[regBInfo['feature_type'] == 'feature_type=Promoter'].shape[0] > 0
		enhancerPresentB = regBInfo[regBInfo['feature_type'] == 'feature_type=Enhancer'].shape[0] > 0
        
		# add regulatory label to df as appropriate (EP = both enhancer AND promoter, N = neither)
		if promoterPresentB and enhancerPresentB:
			loops.loc[idx, "regB"] = 'EP'
		elif promoterPresentB:
			loops.loc[idx, "regB"] = 'P'
		elif enhancerPresentB:
			loops.loc[idx, "regB"] = 'E'
		else:
			loops.loc[idx, "regB"] = 'N'

		# note as EP if one contact point contains E/EP and other contains P/EP
		if (enhancerPresentA and promoterPresentB) or (promoterPresentA and enhancerPresentB):
			isEP = True
		else:
			isEP = False

		# note as same TAD if identical non-empty TAD labels
		if loops.loc[idx, 'tadA'] == loops.loc[idx, 'tadB'] and not pd.isnull(loops.loc[idx, 'tadA']):
			sameTAD = True
		else: 
			sameTAD = False
        
		# if either TAD label is empty, mark contact row as TAD-less
		if pd.isnull(loops.loc[idx, 'tadA']) or pd.isnull(loops.loc[idx, 'tadB']):
			isTADless = True
		else:
			isTADless = False

		# add appropriate label to "EP_Status" col in row, N means non-EP
		if isEP and sameTAD:
			loops.loc[idx, "Int_Status"] = "intraEP"
		elif isEP and not sameTAD and not isTADless:
			loops.loc[idx, "Int_Status"] = "interEP"
		elif not isEP and sameTAD:
			loops.loc[idx, "Int_Status"] = "intraN"
		elif not isEP and not sameTAD and not isTADless:
			loops.loc[idx, "Int_Status"] = "interN"
		elif isEP and isTADless:
			loops.loc[idx, "Int_Status"] = "tadlessEP"
		elif not isEP and isTADless:
			loops.loc[idx, "Int_Status"] = "tadlessN"

		# Count contacts for initial report
		if (promoterPresentA and enhancerPresentB) or (enhancerPresentA and promoterPresentB):
			numEP += 1
		if enhancerPresentA and enhancerPresentB:
			numEE += 1
		if promoterPresentA and promoterPresentB:
			numPP += 1
		numLoops += 1

		# Provide initial report about regulatory distribution of loops
		# print "Percent EE: ", float(numEE)/(numLoops)
		# print "Percent PP: ", float(numPP)/(numLoops)
		# print "Percent EP: ", float(numEP)/(numLoops)
		# print "Rows completed: ", numLoops

	# return original loops df with cols (regA and regB) added
	print "Percent EP: ", float(numEP)/(numLoops)
	return loops
	

if __name__ == "__main__":
	# Use CH12 TADs BED file for now
	tadFile = '/data2/josh/expCH12/CH12_lieberman_intra_5kb_domains.bed'
	tads = pd.read_table(tadFile)
	tads.columns = ["chr", "startbp", "endbp", "name", "score", "dir", "startbpB", "endbpB", 'rgb']

	loopFile = '/data2/josh/expCH12/lieberman_loops_CH12.txt'
	loops = pd.read_csv(loopFile, sep="\t", header=None)
	loops.columns = ["chrA", "startbpA", "endbpA", "chrB", "startbpB", "endbpB", "score"]

	tadLabeledLoops = addTadCols(tads, loops)
	tadLabeledLoops.to_csv('/data2/josh/expCH12/tad_labeled_loops.csv', index=True, sep='\t', header=True)

	regFile = '/data2/josh/expCH12/mm9_regulatory_converted.bed'
	regs = pd.read_table(regFile)
	regs.columns = ["chr", "startbp", "endbp", "feature_type", "extra1", "extra2"]

	addEPCols(regs, tadLabeledLoops)
	loops.to_csv('/data2/josh/expCH12/labeled_loops_new.csv', index=True, sep='\t', header=True)
