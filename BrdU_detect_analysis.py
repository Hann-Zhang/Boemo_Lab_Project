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

file_number_list = ['05','06','07']

for file_number in file_number_list:

    f = open('output_barcode{}.detect'.format(file_number),'r')

    sixMerBrdU_list = list()
    ProbBrdU_list = list()

    n = 0
    n1 = 0
    Attempts = 0
    Calls = 0

    readID_list = list()
    chromosome_list = list()
    refStart_list = list()
    refEnd_list = list()
    strand_list = list()
    Calls_list = list()
    Attempts_list = list()
    Prob_Calls_list = list()

    with open('output_barcode{}.detect'.format(file_number), 'r') as fp:
        for count, line in enumerate(fp):
            pass

    Total_count = count 

    for count, line in enumerate(f):

        #ignore the header lines
        if line[0] == '#':
            continue

        #split the line into a list by whitespace
        splitLine = line.rstrip().split()
        
        
        n1 = n

        if line[0] == '>':
                
            n += 1
                
            Calls_store = Calls
            Attempts_store = Attempts
            Calls_list.append(Calls_store)
            Attempts_list.append(Attempts_store)

            readID = splitLine[0][1:]
            chromosome = splitLine[1]
            refStart = int(splitLine[2])
            refEnd = int(splitLine[3])
            strand = splitLine[4]
            
            Calls = 0
            Attempts = 0
            
            readID_list.append(readID)
            chromosome_list.append(chromosome)
            refStart_list.append(refStart)
            refEnd_list.append(refEnd)
            strand_list.append(strand)
                
                    
                    
        else:
            probBrdU = float(splitLine[1])
            Attempts += 1
            
            if probBrdU > 0.5:
                    Calls += 1
        
        
        if n1 != n and n>1:
                
            Prob_Calls = Calls_store/Attempts_store
            Prob_Calls_list.append(Prob_Calls)
        
        
        if count == Total_count:
                
            Prob_Calls = Calls/Attempts
            Calls_list.append(Calls)
            Attempts_list.append(Attempts)
            Prob_Calls_list.append(Prob_Calls)
    
    f.close()
    
    # function to plot the distribution of probability of successful calls
    def Prob_BrdU_plot(Prob_Calls_list):

        # create plot
        fig, ax = plt.subplots()
        # plot histogram for values of signals
        df_divergence_min = min(Prob_Calls_list)
        df_divergence_max = max(Prob_Calls_list)
        ax.set_xlim( df_divergence_min,df_divergence_max)
        #     ax.set_xlim(0,0.1)
        n, bins, patches = ax.hist(x=Prob_Calls_list, bins=20, color='#0504aa', alpha=0.7)
        ax.grid(axis='y', alpha=0.75)

        # set labels and title
        ax.set_xlabel('Probability of BrdU Calls across Reads')
        ax.set_ylabel('Frequency')
        #     ax.set_title()
        fig.tight_layout()

        plt.savefig('BrdU_plot_barcode{}.png'.format(file_number),dpi=300, bbox_inches='tight')
        plt.close()

    Prob_BrdU_plot(Prob_Calls_list)