# File name: labelScript.py
# Author: Joshua Price
# Date created: 04/11/2018
# Date last modified: 04/11/2018
# Python Version: 2.7

# Purpose: Script to generate fully labeled loops .csv and plot relative frequencies

import pybedtools as pbt
import pandas as pd
import genLabeledLoops as gll 
from epInitialAnalysis import determineInterEpRatio as addExplicitLabels

# Locate loop, tad, and reg files
# Loop file should be \t delimited .txt (no header or index)
loopFile = '/data2/josh/expCH12/lieberman_loops_CH12.txt'

# TAD and regulatory build files should be in .bed format with at least first 4 cols
tadFile = '/data2/josh/expCH12/CH12_lieberman_intra_5kb_domains.bed'
regFile = '/data2/josh/expCH12/mm9_regulatory_converted.bed'

# saved .csv will have both header and index
saveFileName = '/data2/josh/expCH12/automated_labeled_loops.csv'

# read loops file and name columns
loops = pd.read_csv(loopFile, sep="\t", header=None)
loops.columns = ["chrA", "startbpA", "endbpA", "chrB", "startbpB", "endbpB", "score"]

# read TAD and regulatory bed files
tads = pbt.BedTool(tadFile)
regs = pbt.BedTool(regFile)

# add (tadA, tadB) columns
tadLabeledLoops = gll.addTadCols(tads, loops)
fullLabeledLoops = gll.addEPCols(regs, tadLabeledLoops)

# add (EP_STATUS) column to df
explicitLoops = addExplicitLabels(fullLabeledLoops)

# save to savseFileName path
explicitLoops.to_csv(saveFileName, index=True, sep='\t', header=True)

# Plot frequencies in this dataset
plotEPFrequencies(explicitLoops)