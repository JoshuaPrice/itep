# File name: epInitialAnalysis.py
# Author: Joshua Price
# Date created: 04/10/2018
# Date last modified: 04/11/2018
# Python Version: 2.7

# Purpose: plots TAD-reg counts using matplotlib
	# plotEPFrequencies:  performs the above
	# autolabel: helper function to add numberical values on bar plot

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plotEPFrequencies(labeledLoopsDf):
	# Purpose: plots TAD-reg counts using matplotlib
	# Inputs: 
		# labeledLoopsDf: pd df (with header) with chrA#, startBpA, endBpA, chrB#, startBpB, endBpB, score, tadA, tadB, regA, regB, EP_Status
	# Outputs:
		# none returned
		# Plots of total inter/intra/non frequencies and EP-differentiated inter/intra/non
	
	# count number of each type of TAD-reg interaction
	numIntraEP = (labeledLoopsDf["Int_Status"]=="intraEP").sum()
	numInterEP = (labeledLoopsDf["Int_Status"]=="interEP").sum()
	numIntraN = (labeledLoopsDf["Int_Status"]=="intraN").sum()
	numInterN = (labeledLoopsDf["Int_Status"]=="interN").sum()
	numTadlessEP = (labeledLoopsDf["Int_Status"]=="tadlessEP").sum()
	numTadlessN = (labeledLoopsDf["Int_Status"]=="tadlessN").sum()
	numTotalValid = numIntraEP + numInterEP + numIntraN + numInterN

	# Determine fractions of TAD-reg interactions
	fracIntraEP = float(numIntraEP)/(numIntraEP+numIntraN)
	fracInterEP = float(numInterEP)/(numInterEP+numInterN)
	fracEPInter = float(numInterEP)/(numIntraEP+numInterEP)
	fracAllInter = float(numInterEP+numInterN)/(numTotalValid)

	# Report fractions of TAD-reg interactions
	print "Fraction of intra-TAD contacts that are EP:", fracIntraEP
	print "Fraction of inter-TAD contacts that are EP:", fracInterEP
	print "Fraction of EP contacts that are inter:", fracEPInter
	print "Fraction of all contacts that are inter:", fracAllInter

	## Plot 1: Bar plot of intra-TAD, inter-TAD, and TAD-less counts
	ind = np.arange(3)
	vals = np.array([numIntraEP+numIntraN,numInterEP+numInterN,numTadlessEP+numTadlessN])
	fig, ax = plt.subplots()
	width = 0.5       # the width of the bars
	rects1 = ax.bar(ind, vals, width, color='k')

	# add some text for labels, title and axes ticks
	ax.set_ylim([0,1800])
	ax.set_xticks(ind)
	ax.set_xticklabels(('Intra-TAD', 'Inter-TAD', 'TAD-less'))
	ax.set_ylabel('Occurrence Count')
	ax.set_title('Spread of all interactions in Rao 2014 CH12 cells')
	autolabel(rects1, ax)


	## Plot 2: Bar plot of intra-TAD, inter-TAD, and TAD-less counts, with E-P vs N distinguished
	ind = np.arange(3)  # the x locations for the groups
	width = 0.35       # the width of the bars
	fig, ax = plt.subplots()

	epVals = np.array([numIntraEP,numInterEP,numTadlessEP])
	epRects = ax.bar(ind, epVals, width, color='k')

	nVals = np.array([numIntraN,numInterN,numTadlessN])
	nRects = ax.bar(ind + width, nVals, width, color='0.5')

	# add some text for labels, title and axes ticks
	ax.set_ylim([0,1100])
	ax.set_ylabel('Occurrence Count')
	ax.set_title('Spread of all interactions in Rao 2014 CH12 cells')
	ax.set_xticks(ind + width / 2)
	ax.set_xticklabels(('Intra-TAD', 'Inter-TAD', 'TAD-less'))
	ax.legend((epRects[0], nRects[0]), ('E-P', 'Not E-P'))
	autolabel(epRects, ax)
	autolabel(nRects, ax)

	# Display plots
	plt.show()


def autolabel(rects, ax):
    """
    Taken from: https://matplotlib.org/examples/api/barchart_demo.html
    Purpose: Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')


if __name__ == "__main__":
	# recall TAD and reg labeled contacts file
	labeledLoopsFile = '/data2/josh/old/expCH12/automated_labeled_loops.csv'
	labeledLoopsDf = pd.read_csv(labeledLoopsFile, sep="\t")
	plotEPFrequencies(labeledLoopsDf)