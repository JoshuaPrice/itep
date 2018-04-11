# File name: epInitialAnalysis.py
# Author: Joshua Price
# Date created: 04/10/2018
# Date last modified: 04/11/2018
# Python Version: 2.7

# Purpose: takes loops df labeled with TADs and reg elements, reports tad-reg information, saves file with explicit tad-reg labels
	# determineInterEPRatio:  performs the above

import pybedtools as pbt
import pandas as pd

def determineInterEpRatio(labeledLoopsDf):
	# Purpose: takes loops df labeled with TADs and reg elements, reports tad-reg information, saves file with explicit tad-reg labels
	# Inputs: 
		# labeledLoopsDf: pd df (with header) with chrA#, startBpA, endBpA, chrB#, startBpB, endBpB, score, tadA, tadB, regA, regB
	# Outputs:
		# labeledLoopsDf: original pd df with one new columns (EP_Status) populated with string corresponding to TAD-reg status

	# Add new columns to loops df for explicit tad-reg combined label
	labeledLoopsDf["EP_Status"] = ""

	for idx, row in labeledLoopsDf.iterrows():
		# note as EP if one contact point contains E/EP and other contains P/EP
		if ((labeledLoopsDf.loc[idx, 'regA'] == 'E' or labeledLoopsDf.loc[idx, 'regA'] == 'EP') and \
			(labeledLoopsDf.loc[idx, 'regB'] == 'P' or labeledLoopsDf.loc[idx, 'regB'] == 'EP')) or \
			((labeledLoopsDf.loc[idx, 'regA'] == 'P' or labeledLoopsDf.loc[idx, 'regA'] == 'EP') and \
			(labeledLoopsDf.loc[idx, 'regB'] == 'E' or labeledLoopsDf.loc[idx, 'regB'] == 'EP')):
			isEP = True
		else:
			isEP = False

		# note as same TAD if identical non-empty TAD labels
		if labeledLoopsDf.loc[idx, 'tadA'] == labeledLoopsDf.loc[idx, 'tadB'] and labeledLoopsDf.loc[idx, 'tadA']:
			sameTAD = True
		else: 
			sameTAD = False

		# if either TAD label is empty, mark contact row as TAD-less
		isTADless = False
		if pd.isnull(labeledLoopsDf.loc[idx, 'tadA']) or pd.isnull(labeledLoopsDf.loc[idx, 'tadB']):
			isTADless = True

		# add appropriate label to "EP_Status" col in row, N means non-EP
		if isEP and sameTAD:
			labeledLoopsDf.loc[idx, "EP_Status"] = "intraEP"
		elif isEP and not sameTAD and not isTADless:
			labeledLoopsDf.loc[idx, "EP_Status"] = "interEP"
		elif not isEP and sameTAD:
			labeledLoopsDf.loc[idx, "EP_Status"] = "intraN"
		elif not isEP and not sameTAD and not isTADless:
			labeledLoopsDf.loc[idx, "EP_Status"] = "interN"
		elif isEP and isTADless:
			labeledLoopsDf.loc[idx, "EP_Status"] = "tadlessEP"
		elif not isEP and isTADless:
			labeledLoopsDf.loc[idx, "EP_Status"] = "tadlessN"
            
    # Count numbers of each type of interaction
	numIntraEP = (labeledLoopsDf["EP_Status"]=="intraEP").sum()
	numInterEP = (labeledLoopsDf["EP_Status"]=="interEP").sum()
	numIntraN = (labeledLoopsDf["EP_Status"]=="intraN").sum()
	numInterN = (labeledLoopsDf["EP_Status"]=="interN").sum()
	numTotal = numIntraEP + numInterEP + numIntraN + numInterN

	# Report ratios of interactions
	print "Fraction of intra-TAD contacts that are EP:", float(numIntraEP)/(numIntraEP+numIntraN)
	print "Fraction of inter-TAD contacts that are EP:", float(numInterEP)/(numInterEP+numInterN)
	print "Fraction of EP contacts that are inter:", float(numInterEP)/(numIntraEP+numInterEP)
	print "Fraction of all contacts that are inter:", float(numInterEP+numInterN)/(numTotal)

	# Return original df with new column providing explicit TAD-reg label
	return labeledLoopsDf

if __name__ == "__main__":
	labeledLoopsFile = '/data2/josh/expCH12/labeled_loops.csv'
	labeledLoops = pd.read_csv(labeledLoopsFile, sep="\t")

	explicitLoops = determineInterEpRatio(labeledLoops)
	explicitLoops.to_csv('/data2/josh/expCH12/labeled_loops_extra.csv', index=True, sep='\t', header=True)
