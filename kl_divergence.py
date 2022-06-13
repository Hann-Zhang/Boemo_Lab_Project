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

kl_barcode_list = list()
file_number_list = ['03','04','05','06','07','08','09','10','11','12']
df_kl_divergence = pd.DataFrame()


for file_number in file_number_list:

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

    f1.close()
    

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
    # which will be called in function later
    exclude_items_sixmer = list()

    for i in exclude_list:
        exclude_items_sixmer.append(sixMer_list[i])
    
    for j in exclude_items_sixmer:
        sixMer_list.remove(j)

    # import model distribution
    f2 = open("6mer_template.model",'r')

    #creates container 
    sixMer_model_list = list()
    mu_model_list = list()
    std_model_list = list()

    lines = f2.readlines()[1:]

    for line in lines:

        #ignore the header lines
        if line[0] == 'kmer':
                continue

        #split the line into a list by whitespace
        splitLine = line.rstrip().split()

        sixMer_model = splitLine[0]
        mu_model = float(splitLine[1])
        std_model = float(splitLine[2])
                
        sixMer_model_list.append(sixMer_model)
        mu_model_list.append(mu_model)
        std_model_list.append(std_model)
                
        #add these values to a container or do some processing here

    f2.close()

    #creates dictonary for sixmer and model parameters 
    d_model ={}
    for i in range(len(sixMer_model_list)):
        d_model[sixMer_model_list[i]] = [mu_model_list[i],std_model_list[i]]
    
    # function to calculate KL-divergence of p from q given parameters of gaussian distributions
    def KL_divergence(mu_1,sigma_1,mu_2,sigma_2):
        kl = np.log(sigma_2/sigma_1) + (sigma_1**2 +(mu_1-mu_2)**2)/(2*sigma_2**2) -1/2
        return kl

    # function to calculate KL-divergence for each sixmer 
    def KL_divergence_calculate(k):
        
        #gaussian mixture parameter
        mu_1 = d_GMM[k][0]
        sigma_1 = d_GMM[k][1]
        mu_2 = d_GMM[k][2]
        sigma_2 = d_GMM[k][3]
        
        #model parameter
        mu_model = d_model[k][0]
        sigma_model = d_model[k][1]
        
        # calculating kl-divergence 
        kl1 = KL_divergence(mu_1,sigma_1,mu_model,sigma_model)
        kl2 = KL_divergence(mu_2,sigma_2,mu_model,sigma_model)
        
        #assuming the distribution with mean farthest away from the model mean was the analogue 6mer
        #take bigger kl divergence as the analogue 6mer divergence - defined as kl2 (second distribution)
        #take smaller kl divergence as the distribution match to the ONT model - defined as kl1 (first distribution)
        kl_1 = min(kl1,kl2)
        kl_2 = max(kl1,kl2)
        
        return [kl_1,kl_2]
    

    # calculate KL-divergence for each 6mer in sample dataset
    kl_list = list()

    # only for sixmers containing analogue base (sixMer_list)  
    for k in sixMer_list:
        kl = KL_divergence_calculate(k)
        kl_list.append(kl)

    kl_1_list=[]
    kl_2_list=[]
    for i in range(len(kl_list)):
        kl_1_list.append(kl_list[i][0])
        kl_2_list.append(kl_list[i][1])

    df_kl_divergence = pd.DataFrame()
    df_kl_divergence['barcode_{}_01'.format(file_number)]=kl_1_list
    df_kl_divergence['barcode_{}_02'.format(file_number)]=kl_2_list

    df_kl_divergence.to_csv('kl_divergence_barcode_{}.csv'.format(file_number))  
 
