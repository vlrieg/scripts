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
    samples = [] #list of k values
    binom_output = [] #list of the different probabilities of k successes calculated by sample_binom()

    for i in range(sim_num):
        output_tuple = sample_binom(n_num, prob_num)
        samples.append(output_tuple[0])
        binom_output.append(output_tuple[1])
    
    return(samples, binom_output)

#create a list of different probabilities
#from 0.05 to 0.95 with an interval of 1 (0.05...0.15...0.25...etc)
prob_list = np.linspace(0.05, 0.95, 10)
#short_list = [0.05]

for number in prob_list:
#for number in short_list:
    #set input values
    trials = 10000
    n_val = 30
    prob_val = number
    
    #run the function
    sims = many_k_vals(trials, n_val, prob_val)

    #store function output in new variables for simplicity  
    k_val_list = sims[0]
    binom_output_list = sims[1]

    #calculate the sum and average
    sample_sum = sum(k_val_list)
    sample_ave = sample_sum / len(k_val_list)
    print(f"sum of all ks = {sample_sum}, ave k value = {sample_ave}, and probability = {prob_val}")   
    #print(f"sum of all ks = {sample_sum}, ave k value = {sample_ave}, and binom function output = {binom_output_list}")   

    plt.scatter(k_val_list, binom_output_list, label=f"{prob_val}")
    #plt.plot(k_val_list, binom_output_list, label=f"{prob_val}")


plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=9,
           ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel('k value')
plt.ylabel('Outcome Probability')
plt.title(f"The most probable k value for {n_val} trials with varying probabilities of success")

plt.savefig(f"MostLikelyKvalue_for{n_val}trials_varyingprobability.pdf", bbox_inches = "tight")



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
