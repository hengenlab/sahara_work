import musclebeachtools_hlab as mbt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from datetime import datetime as dt
import time
import glob, os
import sahara_work as s
import pickle
import pandas as pd

def all_cv(spikes, step, start_age, end_age):
    edges = np.arange(start_age, end_age, step)
    cvs = np.repeat(np.nan, len(edges)-1)
    for i in np.arange(0,np.size(edges)-1):
        tmp = cv_isi(edges[i], edges[i+1], spikes)
        cvs[i] = tmp
    return cvs, edges

def cv_isi(start, stop, times):
    spks = times[np.where(np.logical_and(times>=start, times<stop))]
    isi = np.diff(spks)
    cv = np.std(isi)/np.mean(isi)
    return cv

bigtic = time.time()
geno = ['te4']

paths = s.get_paths(scorer = 'xgb', geno = geno)

# Subset of the full dataset for testing (comment out when not needed)
#paths = np.load('/media/HlabShare/AD_paper/fr_subset_paths.npy')

# for preallocation, figure out a reasonable number larger than the total number of cells. 100 in a block would be really high yield, the average will definitely be lower than that. 
num_paths = len(paths)
cellnum = num_paths * 100
cvstep = 300

maxt = 12*3600 # seconds in a clustering output
tic = time.time()
print('---- allocating array space')
fr_array = np.zeros(maxt*cellnum, dtype=np.int16)
fr_array[fr_array==0]=-1

cellid_array = np.zeros(maxt*cellnum, dtype=np.int16)
cellid_array[cellid_array==0]=-1

time_array = np.zeros(maxt*cellnum, dtype=np.int32)
time_array[time_array==0]=-1

animal_array = np.zeros(maxt*cellnum, dtype=np.int8)
animal_array[animal_array==0]=-1

cellcount_array = np.zeros(maxt*cellnum, dtype=np.int32)
cellcount_array[cellcount_array==0]=-1

cellqual_array = np.zeros(maxt*cellnum, dtype=np.int8)
cellqual_array[cellqual_array==0] = -1

cv_array = np.zeros(int(maxt/cvstep)*cellnum)
cv_array[cv_array == 0] = -1

cv_cellid_array = np.zeros(int(maxt/cvstep)*cellnum, dtype = np.int16)
cv_cellid_array[cv_cellid_array == 0] = -1

cv_time_array = np.zeros(int(maxt/cvstep)*cellnum, dtype = np.int32)
cv_time_array[cv_time_array == 0] = -1

cv_animal_array = np.zeros(int(maxt/cvstep)*cellnum, dtype = np.int8)
cv_animal_array[cv_animal_array == 0] = -1

cv_cellcount_array = np.zeros(int(maxt/cvstep)*cellnum, dtype = np.int32)
cv_cellcount_array[cv_cellcount_array == 0] = -1

cv_cellqual_array = np.zeros(int(maxt/cvstep)*cellnum, dtype = np.int8)
cv_cellqual_array[cv_cellqual_array == 0] = -1

toc = time.time()
print(f'time to allocate space: {(toc-tic)/60} min')

tic = time.time()
print("-----collecting fr data")
cellcount = 0
tcount = 0
cvcount = 0
for idx, p in enumerate(paths):
    print(f'{idx} of {len(paths)}')
    animal, _, _, _ = s.get_info_from_path(p)
    cells = np.load(p, allow_pickle=True)
    age_sec = s.get_age_sec(start_time = cells[0].rstart_time, birthday = s.get_birthday(animal))
    for cell in cells:
        if cell.quality < 4:
            if cell.quality == 0:
                print('found a qual 0 cell')
            fr_edges = np.arange(cell.start_time + age_sec, cell.end_time + age_sec, 1)
            vals, bins = np.histogram(cell.spike_time_sec+age_sec, fr_edges)
            fr_array[tcount:tcount+len(vals) ] = vals
            cellid_array[tcount:tcount+len(vals) ] = np.repeat(cell.clust_idx,len(vals))
            time_array[tcount:tcount+len(vals) ] = bins[0:-1]
            animal_array[tcount:tcount+len(vals) ] = np.repeat(s.encode_animal(animal), len(vals))
            cellqual_array[tcount:tcount+len(vals) ] = np.repeat(cell.quality, len(vals))
            cellcount_array[tcount:tcount+len(vals) ]= np.repeat(cellcount, len(vals))

            cv, edges = all_cv(cell.spike_time_sec+age_sec, cvstep, cell.start_time+age_sec, cell.end_time+age_sec)
            cv_array[cvcount:cvcount+len(cv)] = cv
            cv_cellid_array[cvcount:cvcount+len(cv)] = np.repeat(cell.clust_idx,len(cv))
            cv_time_array[cvcount:cvcount+len(cv)] = edges[0:-1]
            cv_animal_array[cvcount:cvcount+len(cv)] = np.repeat(s.encode_animal(animal), len(cv))
            cv_cellqual_array[cvcount:cvcount+len(cv)] = np.repeat(cell.quality, len(cv))
            cv_cellcount_array[cvcount:cvcount+len(cv)]= np.repeat(cellcount, len(cv))

            cvcount = cvcount + len(cv)
            tcount = tcount+len(vals)
            cellcount+=1
toc = time.time()
print(f'time to collect {len(paths)} paths worth of data ({cellcount} cells): {(toc-tic)/60} min')
tic = time.time()
print('------deleting nans')
fr_array = np.delete(fr_array, np.where(fr_array == -1))
cellid_array = np.delete(cellid_array, np.where(cellid_array == -1))
time_array = np.delete(time_array, np.where(time_array == -1))
animal_array = np.delete(animal_array, np.where(animal_array == -1))
cellqual_array = np.delete(cellqual_array, np.where(cellqual_array == -1))
cellcount_array = np.delete(cellcount_array, np.where(cellcount_array == -1))

cv_array = np.delete(cv_array, np.where(cv_array == -1))
cv_cellid_array = np.delete(cv_cellid_array, np.where(cv_cellid_array == -1))
cv_time_array = np.delete(cv_time_array, np.where(cv_time_array == -1))
cv_animal_array = np.delete(cv_animal_array, np.where(cv_animal_array == -1))
cv_cellqual_array = np.delete(cv_cellqual_array, np.where(cv_cellqual_array == -1))
cv_cellcount_array = np.delete(cv_cellcount_array, np.where(cv_cellcount_array == -1))

toc = time.time()
print(f'time to delete nans: {(toc-tic)/60} min')

np.save(f'/media/HlabShare/AD_paper/FR_testing/{geno[0]}_fr_array.npy', fr_array.astype(np.int16))
np.save(f'/media/HlabShare/AD_paper/FR_testing/{geno[0]}_cellid_array.npy', cellid_array.astype(np.int16))
np.save(f'/media/HlabShare/AD_paper/FR_testing/{geno[0]}_time_array.npy', time_array.astype(np.int32))
np.save(f'/media/HlabShare/AD_paper/FR_testing/{geno[0]}_animal_array.npy', animal_array.astype(np.int8))
np.save(f'/media/HlabShare/AD_paper/FR_testing/{geno[0]}_cellqual_array.npy', cellqual_array.astype(np.int8))
np.save(f'/media/HlabShare/AD_paper/FR_testing/{geno[0]}_cellcount_array.npy', cellcount_array.astype(np.int32))

np.save(f'/media/HlabShare/AD_paper/FR_testing/DATA/{geno[0]}_cv_array.npy', cv_array)
np.save(f'/media/HlabShare/AD_paper/FR_testing/DATA/{geno[0]}_cv_cellid_array.npy', cv_cellid_array.astype(np.int16))
np.save(f'/media/HlabShare/AD_paper/FR_testing/DATA/{geno[0]}_cv_time_array.npy', cv_time_array.astype(np.int32))
np.save(f'/media/HlabShare/AD_paper/FR_testing/DATA/{geno[0]}_cv_animal_array.npy', cv_animal_array.astype(np.int8))
np.save(f'/media/HlabShare/AD_paper/FR_testing/DATA/{geno[0]}_cv_cellqual_array.npy', cv_cellqual_array.astype(np.int8))
np.save(f'/media/HlabShare/AD_paper/FR_testing/DATA/{geno[0]}_cv_cellcount_array.npy', cv_cellcount_array.astype(np.int32))

bigtoc = time.time()
print(f'TOTAL TIME: {(bigtoc-bigtic)/60} min')