import musclebeachtools_hlab as mbt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from datetime import datetime
import time
import glob, os
import sahara_work as sw
import pickle
import pandas as pd

def cv_isi(start, stop, times):
    spks = times[np.where(np.logical_and(times>=start, times<stop))]
    isi = np.diff(spks)
    cv = np.std(isi)/np.mean(isi)
    return cv

# what's the time window (in seconds) for calculating C.V.?
step = 300

# Load the neuron objects:
wtdat = np.load('wt_compiled.npy', allow_pickle=True)
te4dat = np.load('old_te4_compiled.npy', allow_pickle=True)

cv = []
cellcount = 0
for neuron in wtdat:

    edges = np.arange(neuron.start_time + neuron.age_sec, neuron.end_time + neuron.age_sec, step )

    print(f'working on cell {cellcount}')
    nrn = neuron.spike_time_sec+neuron.age_sec

    for i in np.arange(0,np.size(edges)-1):

        tmp = cv_isi(edges[i], edges[i+1],nrn)

        cv = np.append(cv, [cellcount, edges[i], tmp])

    cellcount+=1

# reshape CV into a n observations by 3 array.
# column 1 is cell count (then you can look at an individual cell/curate)
# column 2 is the bin edge (time in seconds from birth) for the measurement in that row
# column 3 is the C.V. for the cell and bin in a given row
cv = cv.reshape(-1,3)


# ----------------------------------------------------------------
# preallocate to speed things up, convert zeros to NaN so nothing's ambigious.
# Note: this reduces time to <10% of time without preallocation.
maxtimes = np.zeros(len(te4dat))
mcount = 0
for neuron in te4dat:
    maxtimes[mcount] = neuron.end_time
    mcount+=1
maxt = np.int(np.ceil(np.max(maxtimes)))
fr_array = np.zeros(maxt*len(te4dat))
fr_array[fr_array==0]=np.nan
cellid_array = np.zeros(maxt*len(te4dat))
cellid_array[cellid_array==0]=np.nan
time_array = np.zeros(maxt*len(te4dat))
time_array[time_array==0]=np.nan
animal_array = np.zeros(maxt*len(te4dat))
animal_array[animal_array==0]=np.nan 


a = time.time()

cellcount = 0
tcount = 0
for neuron in te4dat:

    fr_edges = np.arange(neuron.start_time + neuron.age_sec, neuron.end_time + neuron.age_sec, 1 )

    tmp = np.histogram(neuron.spike_time_sec+neuron.age_sec, fr_edges)

    fr_array[tcount:tcount+len(tmp[0]) ] = tmp[0]
    cellid_array[tcount:tcount+len(tmp[0]) ] = np.repeat(cellcount,len(tmp[0]))
    time_array[tcount:tcount+len(tmp[0]) ] = tmp[1][0:-1]
    animal_array[tcount:tcount+len(tmp[0]) ]= np.repeat(sw.encode_animal(te4dat[0].animal), len(tmp[0]) )

    tcount = tcount+len(tmp[0])
    cellcount+=1

b = time.time()
print(f'Elapsed time for FR repackaging is {b-a} seconds.')

#remove the NaNs at the end of the arrays:
fr_array = np.delete(fr_array, np.where(np.isnan(fr_array)))
cellid_array = np.delete(cellid_array, np.where(np.isnan(cellid_array)))
time_array = np.delete(time_array, np.where(np.isnan(time_array)))
animal_array = np.delete(animal_array, np.where(np.isnan(animal_array)))

np.save(fn + 'te4firing.npy', fr_array.astype(np.int16))
np.save(fn + 'te4cellid.npy', cellid_array.astype(np.in32))
np.save(fn + 'te4time.npy', time_array.astype(np.int32))
np.save(fn + 'te4animal.npy', animal_array.astype(np.int8))

c = time.time()
print(f'All told, that took {c-a} seconds to run and save.')