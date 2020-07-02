import numpy as np
from sahara_work import Criticality_new as cr
from copy import deepcopy as cdc
import powerlaw as plaw 
from collections import Counter

def initial_exclude(burst, t, bound=4, freq=20):
    '''
        this function runs on the logic outlined in beggs 2016 for doubly-truncated power-laws

        it discards any avalanch that has a size OR durration of less than 4.
        then it discards any avalanch with a size OR durration that occurs less than 20 times

        for example, if a size of 12 appears only 7 times in the whole dataset, any avalanch with a size of 12 
            is discarded even if the corresponding durration appears more than 20 times



        params:
            burst distribution
            time distribution
            bound: boundary for the minimum. we're using 4 as that's what beggs used
            freq: boundary for the frequency of appearance, we're using 20 as that's what beggs used

        returns:
            newly pruned dataset to then do further truncating on for power-law fits
            burst, t
    '''
    # first get locations of burst and durrations less than 4
    locs_boundary = np.where(np.logical_or(burst < 4, t < 4))[0]

    # get counts for each occurance 
    counts_bursts = Counter(burst)
    counts_t = Counter(t)

    #range(len(burst)) is an array of all the indicies in burst and t since they're the same size
    #checks to see if each of the values in burst or t occurs less than 20 times and records their indicies if they do
    locs_to_rem = [i for i in range(len(burst)) if counts_burst[burst[i]] < 20 or counts_t[t[i]] < 20]

    #combines two location arrays
    final_locs_to_remove = np.union1d(locs_boundary, locs_to_rem)

    #turns the indicies into true/false values 
    indicies = np.isin(range(len(burst)), final_locs_to_remove)

    #takes everything other than the idicies specified above
    b = burst[~indicies]
    t = t[~indicies]

    return b, t










