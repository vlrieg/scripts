#!/usr/bin/env python3

import math
import sys
import random

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np

##############################################################################
# binomial distribution function
def binom(n, k, p):
    # calculate number of arrangements
    order = math.factorial(n) / (math.factorial(k) * math.factorial(n - k))
    #print(order)

    # calculate probability of success
    prob_success = (p ** k) * ((1 - p) ** (n - k))
    #print(prob_success)
    
    # calculate the outcome
    prob_n_choose_k = order * prob_success
    #print(prob_n_choose_k)
    return(prob_n_choose_k)


#n_input = int(sys.argv[1]) # number of trials
#k_input = int(sys.argv[2]) # number of successes
#p_input = float(sys.argv[3]) # probability of success in each trial    
##############################################################################

# random k value
def sample_binom(n, p):
    u = random.uniform(0,1)                                
    sum=0.0

    assert p <= 1, "P is not less or equal to 1!!!!"


    for kvalue in range(n+1):
        pr = binom(n,kvalue,p)
        #print("k={} pr={}".format(kvalue,pr))
        sum += pr
        if sum >= u:
            #print(pr) #this prints the pr output value successfully
            return(kvalue, pr)

    raise "Oh no I didn't stop!"

##############################################################################
# create a function that selects a random k value many different times 
# using the sample_binom() function
# The # in range(#) is the number of times a k value is "simulated"

def many_k_vals(sim_num, n_num, prob_num):
    tuple_list = []

    for i in range(sim_num):
        output_tuple = sample_binom(n_num, prob_num)
        tuple_list.append(output_tuple)
    
    return(tuple_list)

#create a list of different probabilities
#from 0.05 to 0.95 with as many intervals as designated by num_intervals variable
num_intervals = 10
prob_list = np.linspace(0.05, 0.95, num_intervals)
#short_list = [0.05]

for number in prob_list:
#for number in short_list:
    #set input values
    trials = 10000
    n_val = 8
    prob_val = number
    
    #run the function
    sims = many_k_vals(trials, n_val, prob_val) #sims now contains a list of tuples, each tuple = (x, y)

    sorted_sims = sorted(sims, key=lambda tuple: tuple[0]) #sort on the first index of each tuple (x value)

    #store sorted many_k_vals() output in new lists: one for x values and one for y values
    k_val_list = [item[0] for item in sorted_sims] #x values
    binom_output_list = [item[1] for item in sorted_sims] #y values

    #calculate the sum and average
    sample_sum = sum(k_val_list)
    sample_ave = sample_sum / len(k_val_list)
    #print(f"sum of all ks = {sample_sum}, ave k value = {sample_ave}, and probability = {prob_val}")   
    #print(f"sum of all ks = {sample_sum}, ave k value = {sample_ave}, and binom function output = {binom_output_list}")   

    #plt.scatter(k_val_list, binom_output_list, label="{0:.4f}".format(prob_val))
    #plt.plot(k_val_list, binom_output_list, label=f"{prob_val}") # will print long list of decimal places in legend
    plt.plot(k_val_list, binom_output_list, label="{0:.4f}".format(prob_val)) #limit to 4 decimal places using formatting

    

    
plt.legend(loc="upper left", bbox_to_anchor=(1,1))
plt.xlabel('k value')
plt.ylabel('Outcome Probability')
plt.title(f"The most probable k value for {n_val} trials with varying probabilities of success")

plt.savefig(f"MostLikelyKvalue_for{n_val}trials_varyingprobability_{num_intervals}probintervals.pdf", bbox_inches = "tight")

# plot the most likely k value for number_n trials and prob_p probability
#plt.scatter(k_val_list, binom_output_list)

#plt.xlabel('k value')
#plt.ylabel('Outcome Probability')
#plt.title(f"The most probable k value for {n_val} trials with varying probabilities of success")

#plt.savefig(f"MostLikelyKvalue_for{n_val}trials_{prob_val}probability.pdf")
##############################################################################


# # Compute the expected value (aka the limit of the means) using the sum pk*k formula
# pk_times_k = [] #empty list

# for kval in samples:
#     pk_times_k.append((kval/sample_sum)*kval) #kval/sample_sum gives the probability of each k value


# sum_kprob_timesk = sum(pk_times_k)
# print(f"The expected k value is {sum_kprob_timesk}.")



##############################################################################
# Hold n and p constant but vary the value of k
# p is the "final" probability of exactly k successes in n independent experiments
# for kvalue in range(n_input+1):
#     p = binom(n_input,kvalue,p_input)
#     print("k={} p={}".format(kvalue,p))
#     sum += p
#     if sum > u and chosen is None:
#         chosen = kvalue
    

#print(u)
#print("chosen = {}".format(chosen))    
#print("sum = {}".format(sum))

#binom(n_input, k_input, p_input)
