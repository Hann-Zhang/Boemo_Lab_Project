#! /usr/bin/python

import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
import os

from scipy.stats import norm
from numpy.random import seed
from numpy.random import normal

seed(1)
# dataset2 barcodes: '03','04','05','06','07','08','09','10','11','12'
file_number_list = ['12']

for file_number in file_number_list:

    #load file from DNAscent align

    f = open("output_barcode{}.align".format(file_number),'r')

    #creates container 
    readID_list = list()
    chromosome_list = list()
    refStart_list = list()
    refEnd_list = list()
    strand_list = list()
    sixMerOnRef_list = list()
    signal_list = list()


    for line in f:

        #ignore the header lines
        if line[0] == '#':
                continue

        #split the line into a list by whitespace
        splitLine = line.rstrip().split()

        if line[0] == '>':

                readID = splitLine[0][1:]
                chromosome = splitLine[1]
                refStart = int(splitLine[2])
                refEnd = int(splitLine[3])
                strand = splitLine[4]
                
                readID_list.append(readID)
                chromosome_list.append(chromosome)
                refStart_list.append(refStart)
                refEnd_list.append(refEnd)
                strand_list.append(strand)
                
        else:
                sixMerOnRef = splitLine[4]
                signal = float(splitLine[2])
                    
                sixMerOnRef_list.append(sixMerOnRef)
                signal_list.append(signal)
                
                #add these values to a container or do some processing here

    f.close()

    # create dataframe for sixmer-signal
    df = pd.DataFrame({
        'sixMer': sixMerOnRef_list,
        'signal': signal_list
    })

    # group signals for the same sixMer
    df1 = df.groupby('sixMer').agg(list)
    # extract values of 95% and 5% quantile 
    df1['maxx'] = df.groupby('sixMer').quantile(0.995)
    df1['minn'] = df.groupby('sixMer').quantile(0.005)

    #load file from DNAscent trainGMM
    f1 = open("output_barcode{}.trainGMM".format(file_number),'r')

    #creates container 
    sixMer_list = list()
    mu_1_list = list()
    std_1_list = list()
    mu_2_list = list()
    std_2_list = list()

    lines = f1.readlines()[1:]

    for line in lines:

        #ignore the header lines
        if line[0] == '6mer':
                continue

        #split the line into a list by whitespace
        splitLine = line.rstrip().split()

        sixMer = splitLine[0]
        mu_1 = float(splitLine[4])
        std_1 = float(splitLine[5])
        mu_2 = float(splitLine[7])
        std_2 = float(splitLine[8])
                
        sixMer_list.append(sixMer)
        mu_1_list.append(mu_1)
        std_1_list.append(std_1)
        mu_2_list.append(mu_2)
        std_2_list.append(std_2)
                
        #add these values to a container or do some processing here

    f.close()

    #creates dictonary for sixmer and GMM parameters 
    d_GMM ={}
    for i in range(len(sixMer_list)):
        d_GMM[sixMer_list[i]] = [mu_1_list[i],std_1_list[i],mu_2_list[i],std_2_list[i]]

    # create list with values needs to be excluded in the analysis
    if file_number in ['03','04','05','06','07','08']:
        # remove sixmers without including sixMer analogue base
        exclude_list = list()
        for n, sixMer in enumerate(sixMer_list):
            if "T" not in sixMer:
                exclude_list.append(n)
        
    if file_number in ['11','12']:
        # remove sixmers without including sixMer analogue base
        exclude_list = list()
        for n, sixMer in enumerate(sixMer_list):
            if "A" not in sixMer:
                exclude_list.append(n)
    
    if file_number == '09':
        # remove sixmers without including sixMer analogue base
        exclude_list = list()
        for n, sixMer in enumerate(sixMer_list):
            if "G" not in sixMer:
                exclude_list.append(n)
    
    if file_number == '10':
        # remove sixmers without including sixMer analogue base
        exclude_list = list()
        for n, sixMer in enumerate(sixMer_list):
            if "C" not in sixMer:
                exclude_list.append(n)

    # create a vector with sixmers containing analogue base 
    # which will be called in later function
    exclude_items_sixmer = list()

    for i in exclude_list:
        exclude_items_sixmer.append(sixMer_list[i])
    
    for j in exclude_items_sixmer:
        sixMer_list.remove(j)


    #function to plot fitted Gaussian distributions and raw data 

    def genPlot(k):

        # there are cases where 6mers in trainGMM file 
        # is not present in alignment output
        if k not in df1.index:
            print(f'{k} not present')
            return None

        # signal values belong to the 'k' sixmer
        signalRow = df1.loc[[k]]['signal']
        signalXlimMin = df1.loc[[k]]['minn'].values[0]
        signalXlimMax = df1.loc[[k]]['maxx'].values[0]
        # quantities of signals for AAAAAG sixmer
        elCount = len(signalRow.values[0])
        # full range of signals 
        fullSize = max(signalRow.values[0]) - min(signalRow.values[0])
        # range for quantiles 5%-95%
        scaledSize = signalXlimMax - signalXlimMin
        # set the number of bins for plotting
        nBinsInScreen = 200
        nBins = int(200 * fullSize / scaledSize)
        
        # create plot
        fig, ax = plt.subplots()
        # plot histogram for values of signals
        ax.set_xlim(signalXlimMin,signalXlimMax)
        n, bins, patches = ax.hist(x=signalRow, bins=nBins, color='#0504aa',
                                    alpha=0.7)
        ax.grid(axis='y', alpha=0.75)
        
        # plot the curves for normal distributions 
        mu_1 = d_GMM[k][0]
        sigma_1 = d_GMM[k][1]
        #y1 = ((1 / (np.sqrt(2 * np.pi) * sigma_1)) * np.exp(-0.5 * (1 / sigma_1 * (bins - mu_1))**2))
        #ax.plot(bins, y1, '--',color='green')
        p1 = norm.pdf(bins, mu_1, sigma_1)
        ax.plot(bins, p1 * elCount/p1.sum(), color='green')
        
        mu_2 = d_GMM[k][2]
        sigma_2 = d_GMM[k][3]
        #y2 = ((1 / (np.sqrt(2 * np.pi) * sigma_2)) * np.exp(-0.5 * (1 / sigma_2 * (bins - mu_2))**2))
        #ax.plot(bins, y2, '--',color='red')
        p2 = norm.pdf(bins, mu_2, sigma_2)
        ax.plot(bins, p2 * elCount/p2.sum(), color='red')
        
        # set labels and title
        ax.set_xlabel('signals')
        ax.set_ylabel('Frequency')
        ax.set_title(f'{k} ')
        fig.tight_layout()
        
        # save the plots into directory
        import os
        script_dir = os.path.dirname('/Users/hanzhang/Desktop/output_01062022_DNAscent/')
        results_dir = os.path.join(script_dir, "plots_barcode{}/".format(file_number))
        file_name = '{}.png'.format(k)
        
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        plt.savefig(results_dir+file_name,dpi=300, bbox_inches='tight')
        
        plt.close()
        
    # plot for n number of sixmers
    for n,k in enumerate(sixMer_list):
        genPlot(k)
        # if n > 50:
        #     break



