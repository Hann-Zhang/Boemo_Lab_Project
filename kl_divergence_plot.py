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

file_number_list = ['03','04','05','06','07','08','09','10','11','12']

for file_number in file_number_list:
    # import saved csv file as dataframe
    df_divergence = pd.read_csv('kl_divergence_barcode_{}.csv'.format(file_number),index_col=0)

    def divergence_plot(df_divergence,file_number,kl_number):
        df_divergence_max = df_divergence['barcode_{}_{}'.format(file_number,kl_number)].quantile(0.95)
        df_divergence_min = df_divergence['barcode_{}_{}'.format(file_number,kl_number)].quantile(0)
        
        # full range of signals 
        fullSize = max(df_divergence['barcode_{}_{}'.format(file_number,kl_number)]) - min(df_divergence['barcode_{}_{}'.format(file_number,kl_number)])
        # range for quantiles 5%-95%
        scaledSize = df_divergence_max - df_divergence_min 
        # set the number of bins for plotting
        nBinsInScreen = 20
        nBins = int(20 * fullSize / scaledSize)
        
        # create plot
        fig, ax = plt.subplots()
        # plot histogram for values of signals
        ax.set_xlim( df_divergence_min,df_divergence_max)
        #ax.set_xlim(0,1)
        n, bins, patches = ax.hist(x=df_divergence['barcode_{}_{}'.format(file_number,kl_number)], bins=nBins, color='#0504aa',
                                    alpha=0.7)
        ax.grid(axis='y', alpha=0.75)
        
        # set labels and title
        ax.set_xlabel('KL divergence')
        ax.set_ylabel('Frequency')
        ax.set_title('barcode_{}_kl_{}'.format(file_number,kl_number))
        fig.tight_layout()
        
        # save the plots into directory 
        import os
        script_dir = os.path.dirname('/Users/hanzhang/Desktop/output_01062022_DNAscent/')
        results_dir = os.path.join(script_dir, "plots_KL_divergence_barcode{}/".format(file_number))
        file_name = 'kl_divergence_barcode_{}_kl_{}.png'.format(file_number,kl_number)
        
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        plt.savefig(results_dir+file_name,dpi=300, bbox_inches='tight')
        
        plt.close()
        

    # generate plots

    for kl_number in ['01','02']:
        divergence_plot(df_divergence,file_number,kl_number)

