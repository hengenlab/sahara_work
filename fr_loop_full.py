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

geno = ['te4']

ps = glob.glob('/media/HlabShare/clayton_sahara_work/clustering/*/*/*/*/co/*scored_xgb*.npy')
paths = []
for p in ps:
    animal, _, _, _ = s.get_info_from_path(p)
    if s.get_genotype(animal) in geno:
        paths.append(p)

num_paths = len(paths)
cellnum = num_paths * 100

maxt = 12*3600 # seconds in a clustering output
tic = time.time()
print('---- allocating array space')
fr_array = np.zeros(maxt*cellnum)
fr_array[fr_array==0]=np.nan
cellid_array = np.zeros(maxt*cellnum)
cellid_array[cellid_array==0]=np.nan
time_array = np.zeros(maxt*cellnum)
time_array[time_array==0]=np.nan
animal_array = np.zeros(maxt*cellnum)
animal_array[animal_array==0]=np.nan 
cellcount_array = np.zeros(maxt*cellnum)
cellcount_array[cellcount_array==0]=np.nan
cellqual_array = np.zeros(maxt*cellnum)
cellqual_array[cellqual_array==0] = np.nan

toc = time.time()
print(f'time to allocate space: {(toc-tic)/60} min')

tic = time.time()
print("-----collecting fr data")
cellcount = 0
tcount = 0
for idx, p in enumerate(paths):
    print(f'{idx} of {len(paths)}')
    animal, _, _, _ = s.get_info_from_path(p)
    cells = np.load(p, allow_pickle=True)
    age_sec = get_age_sec(start_time = cells[0].rstart_time, birthday = s.get_birthday(animal))
    for cell in cells:
        if cell.quality < 4:
            
            fr_edges = np.arange(cell.start_time + age_sec, neuron.end_time + age_sec, 1)
            vals, bins = np.histogram(neuron.spike_time_sec+age_sec, fr_edges)
            fr_array[tcount:tcount+len(vals) ] = vals
            cellid_array[tcount:tcount+len(vals) ] = np.repeat(cell.clust_idx,len(vals))
            time_array[tcount:tcount+len(vals) ] = bins[0:-1]
            animal_array[tcount:tcount+len(vals) ]= np.repeat(sw.encode_animal(animal), len(vals) )
            cellqual_array[tcount:tcount+len(vals) ]= np.repeat(cell.quality, len(vals) )
            cellcount_array[tcount:tcount+len(vals) ]= np.repeat(cellcount, len(vals) )

            tcount = tcount+len(vals)
            cellcount+=1
toc = time.time()
print(f'time to collect {len(paths)} paths worth of data ({cellcount} cells): {(toc-tic)/60} min')
tic = time.time()
print('------deleting nans')
fr_array = np.delete(fr_array, np.where(np.isnan(fr_array)))
cellid_array = np.delete(cellid_array, np.where(np.isnan(cellid_array)))
time_array = np.delete(time_array, np.where(np.isnan(time_array)))
animal_array = np.delete(animal_array, np.where(np.isnan(animal_array)))
cellqual_array = np.delete(cellqual_array, np.where(np.isnan(cellqual_array)))
cellcount_array = np.delete(cellcount_array, np.where(np.isnan(cellcount_array)))
toc = time.time()
print(f'time to delete nans: {(toc-tic)/60} min')

np.save('/media/HlabShare/AD_paper/FR_testing/sample_fr_array.npy', fr_array.astype(np.int16))
np.save('/media/HlabShare/AD_paper/FR_testing/sample_cellid_array.npy', cellid_array.astype(np.int16))
np.save('/media/HlabShare/AD_paper/FR_testing/sample_time_array.npy', time_array.astype(np.int32))
np.save('/media/HlabShare/AD_paper/FR_testing/sample_animal_array.npy', animal_array.astype(np.int8))
np.save('/media/HlabShare/AD_paper/FR_testing/sample_cellqual_array.npy', cellqual_array.astype(np.int8))
np.save('/media/HlabShare/AD_paper/FR_testing/sample_cellcount_array.npy', cellcount_array.astype(np.int32))
