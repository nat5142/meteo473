# -------------MileStone 3A------------- #
# Author: Nicholas Tulli
# Email: nat5142@psu.edu
#--------------------------------------- #

import numpy as np
from scipy import stats
from pylab import *
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import random
import UserFunctions as uf

#################################################################################
###############                  Variable Library                 ###############
#################################################################################
B = 0                                                                           #
i = 1                                                                           #
j = 1                                                                           #
k = 1                                                                           #
biased_array = []                                                               #
biased_number = 0                                                               #
dice_numbers = list([1,2,3,4,5,6])                                              #
difference_ratios = []                                                          #
second_array_counts = []                                                        #
biased_counts = []                                                              #
#################################################################################
#################################################################################

#Prompt user to input number of dice roll trials
N_rolls = int(input('How many rolls?:'))
#Generate random sample of N_rolls; the first sample that is generated is an unweighted sample
number_array,array_counts,i = uf.generate_sample(N_rolls,B,i)
#Increase i by one to signal the next step of the generate_sample user function
i+=1

#Generate histogram of the first random sample of dice rolls
number_array,j = uf.plot_histogram(number_array,N_rolls,biased_array,biased_number,B,j,dice_numbers,difference_ratios)
#Increase j by one to signal the next step of the plot_histogram user function
j+=1

#Calculate statistical variables of the first sample of data
first_sd,first_skew,first_kurtosis = uf.sample_stat(array_counts,k,second_array_counts,biased_counts)
#Increase k by one to signal the next step of the sample_stat user function
k+=1

#Prompt user to input new number of dice roll trials to create second sample
N_rolls = int(input('How many rolls for sample 2?:'))
#Prompt user to input a bias in the second sample of dice rolls to create a third, weighted sample
B = int(input('How much bias for sample 2? (enter number between 20 and 100): '))

#Generate both weighted and unweighted sample of N_rolls; the first sample that is generated is an unweighted sample, second is weighted
number_array,biased_array,second_array_counts,biased_counts,i,biased_number = uf.generate_sample(N_rolls,B,i)

#Calculate statistical variables of the second unweighted data sample
second_sd,second_skew,second_kurtosis = uf.sample_stat(array_counts,k,second_array_counts,biased_counts)
#Increase k by one to signal the next step of the sample_stat user function
k+=1
#Calculate statistical variables of the weighted data sample
biased_sd,biased_skew,biased_kurtosis = uf.sample_stat(array_counts,k,second_array_counts,biased_counts)

#Calculate the ratios of the weighted sample to the unweighted sample
difference_ratios = biased_counts / second_array_counts

#Create two-plot figure of data
#first figure shows the frequency of both weighted sample and the second unweighted sample
#second figure displays the ratio of the two samples
number_array,biased_array,j = uf.plot_histogram(number_array,N_rolls,biased_array,biased_number,B,j,dice_numbers,difference_ratios)


























    
    
    
        
