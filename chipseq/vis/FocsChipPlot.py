#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import cooler
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
# plt.style.use('fivethirtyeight')

def genSbsPlot(chipFile, chipName):
    # Import promoter list from FOCS
    ep = pd.read_csv('/data2/josh/ep/fantom_mm10_E-P.txt', sep='\t', header=None, usecols=[0,1,2,3,4,5,6,7])
    ep.columns = ['p_chr', 'p_start', 'p_end', 'p_dir', 'contact_id','e_chr','e_start','e_end']

    # Define centers of E and P and distance
    ep['p_center'] = (ep['p_start']+ep['p_end'])/2
    ep['e_center'] = (ep['e_start']+ep['e_end'])/2
    ep['distance'] = abs(ep['p_center'] - ep['e_center'])
    ep_long = ep[ep['distance'] > 5000]


    # Add chromatin labels

    # Import chip-seq peaks - mm10
    chip = pd.read_csv(chipFile, delimiter='\t', header=None)
    chip.columns = ['chr', 'start', 'end', 'id', 'strength']

    # Label promoters as ChIP-seq overlapping or not (new column)
    promoters_chip_list = [None] * ep.shape[0]
    for idx, row in ep.iterrows():
        chrom = row.p_chr[3:]
        pos = row.p_center
        chip_matches = chip[(chip['chr'] == chrom) & (chip['start'] < pos) & (chip['end'] > pos)]
        if chip_matches.empty:
            promoters_chip_list[idx] = False
        else:
            promoters_chip_list[idx] = True

    ep['p_chip'] = ''
    ep['p_chip'] = promoters_chip_list


    # Label enhancers as ChIP-seq overlapping or not (new column)
    enhancers_chip_list = [None] * ep.shape[0]
    for idx, row in ep.iterrows():
        chrom = row.e_chr[3:]
        pos = row.e_center
        chip_matches = chip[(chip['chr'] == chrom) & (chip['start'] < pos) & (chip['end'] > pos)]
        if chip_matches.empty:
            enhancers_chip_list[idx] = False
        else:
            enhancers_chip_list[idx] = True

    ep['e_chip'] = ''
    ep['e_chip'] = enhancers_chip_list


    # Import and define chromosomes df - mm10
    c = cooler.Cooler('/data2/josh/stan/merge_res200.cool')
    chrs = c.chroms()[:]
    chrs['name'] = chrs['name'].str.replace('chr', '')
    chrs = chrs.set_index('name')
    chrs.index.names = ['chrom']

    # Import TAD list
    tads = pd.read_csv('/data2/josh/tads/TAD_mm10.csv', delimiter='\t')
    tads['chrom'] = tads['chrom'].str.replace('chr','')
    tads['size'] = tads['end'] - tads['start']
    chr_coverage = tads.groupby('chrom').sum()['size'].to_frame()
    chr_coverage = chr_coverage.join(chrs)
    chr_coverage['ratio'] = chr_coverage['size']/chr_coverage['length']
    chr_coverage.sort_values('ratio')


    # tads1['chrom'] = tads1['chrom'].str.replace('chr','')
    # tads1.head()

    # tads1['size'] = tads1['end'] - tads1['start']
    # tads1.groupby('chrom').sum()


    # Add TAD label for promoters

    promoters_tad_list = [None] * ep.shape[0]
    for idx, row in ep.iterrows():
        chrom = row.p_chr
        pos = row.p_center
        tad_matches = tads[(tads['chrom'] == chrom) & (tads['start'] < pos) & (tads['end'] > pos)]
        if tad_matches.empty:
            promoters_tad_list[idx] = -1
        else:
            promoters_tad_list[idx] = tad_matches.index[0]

    ep['p_tad'] = ''
    ep['p_tad'] = promoters_tad_list


    # Add TAD label for enhancers

    enhancers_tad_list = [None] * ep.shape[0]
    for idx, row in ep.iterrows():
        chrom = row.e_chr
        pos = row.e_center
        tad_matches = tads[(tads['chrom'] == chrom) & (tads['start'] < pos) & (tads['end'] > pos)]
        if tad_matches.empty:
            enhancers_tad_list[idx] = -1
        else:
            enhancers_tad_list[idx] = tad_matches.index[0]

    ep['e_tad'] = ''
    ep['e_tad'] = enhancers_tad_list


    # Save TAD- and ChIP-labeled E-P labels to bed file
    ep.to_csv('labeled_ep.bed',sep='\t',index=True, header=True)

    # If file already saved, read in instead
    # ep=pd.read_csv('labeled_ep.bed',sep='\t')


    # Count numbers in each category: E-enhriched, P-enriched, EP-enriched for inter- and intra-TAD contacts
    # Only consider pairs where both E and P are in a labeled TAD
    ep_in_tads = ep[(ep['p_tad']!=-1) & (ep['e_tad']!=-1)]


    # define helper function tad_status
    def tad_status(p_tad, e_tad):
        if (p_tad == e_tad):
            return 'Intra-TAD'
        else:
            return 'Inter-TAD'
        
    # Create list for EP_status the put it into the df 
    tad_status_list = [tad_status(row['p_tad'], row['e_tad']) for idx, row in ep_in_tads.iterrows()]
    ep_in_tads['tad_status'] = tad_status_list


    # define helper function ep_status
    def ep_status(e_chip, p_chip):
        if ((e_chip == True) and (p_chip == True)):
            return 'EP'
        elif ((e_chip == False) and (p_chip == True)):
            return 'P'
        elif ((e_chip == True) and (p_chip == False)):
            return 'E'
        else:
            return 'N'

    # Create list for EP_status the put it into the df 
    ep_status_list = [ep_status(row['e_chip'], row['p_chip']) for idx, row in ep_in_tads.iterrows()]
    ep_in_tads['ep_status'] = ep_status_list


    # Separate df into two dfs (one inter, one intra) for plot generation
    intra_tad_ep = ep_in_tads[ep_in_tads['tad_status']=='Intra-TAD']
    inter_tad_ep = ep_in_tads[ep_in_tads['tad_status']=='Inter-TAD']


    # Generate sbs (side-by-side) plot
    plt.figure(figsize=[10,8])
    plt.subplot(221)

    ax1 = sns.countplot(x="ep_status", data=intra_tad_ep)
    plt.xlabel('Type of Interaction')
    plt.ylabel('Number of contacts')
    plt.title('Intra-TAD EPIs')
    L=plt.legend()

    plt.subplot(222)
    ax2 = sns.countplot(x="ep_status", data=inter_tad_ep)
    plt.xlabel('Type of Interaction')
    plt.ylabel('Number of contacts')
    plt.title('Inter-TAD EPIs')
    L=plt.legend()
    plt.suptitle(chipName + ' Enrichment in Intra- and Inter-TAD EPIs',size=20)

    # Save figure
    plt.savefig(chipName + '_sbs.png')

    # Calculate exact proportions

    # Count and calculate proporiton of intra-TAD E-enriched, P-enriched, EP-enriched
    p_enriched_intra = intra_tad[(intra_tad['p_chip']==True) & (intra_tad['e_chip']==False)].shape[0]
    e_enriched_intra = intra_tad[(intra_tad['p_chip']==False) & (intra_tad['e_chip']==True)].shape[0]
    ep_enriched_intra = intra_tad[(intra_tad['p_chip']==True) & (intra_tad['e_chip']==True)].shape[0]

    p_enriched_intra_proportion = float(p_enriched_intra)/intra_count
    print('P-enriched propotion: ', p_enriched_intra_proportion)

    e_enriched_intra_proportion = float(e_enriched_intra)/intra_count
    print('E-enriched proportion: ', e_enriched_intra_proportion)

    ep_enriched_intra_proportion = float(ep_enriched_intra)/intra_count
    print('EP-enriched count: ', ep_enriched_intra_proportion)


    # Count inter-TAD EP contacts
    inter_tad = ep_in_tads[ep_in_tads['p_tad'] != ep_in_tads['e_tad']]


    # Count and calculate proporiton of inter-TAD E-enriched, P-enriched, EP-enriched
    p_enriched_inter = inter_tad[(inter_tad['p_chip']==True) & (inter_tad['e_chip']==False)].shape[0]
    e_enriched_inter = inter_tad[(inter_tad['p_chip']==False) & (inter_tad['e_chip']==True)].shape[0]
    ep_enriched_inter = inter_tad[(inter_tad['p_chip']==True) & (inter_tad['e_chip']==True)].shape[0]

    p_enriched_inter_proportion = float(p_enriched_inter)/inter_count
    print('P-enriched propotion: ', p_enriched_inter_proportion)

    e_enriched_inter_proportion = float(e_enriched_inter)/inter_count
    print('E-enriched proportion: ', e_enriched_inter_proportion)

    ep_enriched_inter_proportion = float(ep_enriched_inter)/inter_count
    print('EP-enriched count: ', ep_enriched_inter_proportion)


if __name__ == "__main__":
    # recall TAD and reg labeled contacts file
    genSbsPlot('/data2/josh/chipseq/H3K4ME3/macs14_peaks.bed', 'H3K4me3')