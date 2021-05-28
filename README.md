## Crit_utils
These are functions to make running large amounts of data through criticality easier and more efficient in addition to helper functions to deal with data.

## Crit_class
Class structure for running and storing criticality results

## Crit_hlab
Inherits Crit_class, with some added functionality specific to our data structure and functions that use these objects (lilo_and_stitch)

---
## lilo_and_stitch
Outer function that creates and saves Crit objects. Most complicated of the functions in crit_utils. 

This bad boy takes in a list of paths, in addition to a dictionary of parameters, and runs criticality for all the paths and saves criticality objects for each run.

**Parameters**:


*Paths*: a list of paths to each of the scrubbed neurons you want to run criticality on. Make sure the path is **complete**, and goes all the way to a .npy file. 
ex: paths = ['/media/Hlabshare/Clustering_Data/CAF00019/caf19_03262020/0_12/*neurons_group0.npy', ...]


*Params*:
A dictionary with the necessary parameters, there are a many necessary inputs so the function can run. The params dictionary needs look *at least* like this:

```python
params = {
    'flag': 2,  # 1 is DCC, 2 is p_val and DCC
    'ava_binsz': 0.04,  # in seconds
    'hour_bins': 4,  # durration of block to look at. One DCC value will be returned for every 4 hours, e.g.
    'perc': 0.35, # threshold percentage for determining avalanche boarders
    'bm':None, # upper limit of the size min cut off
    'tm':None, # upper limit of the time min cut off
    'nfactor_bm': 0, # where to start the search for a size min
    'nfactor_tm': 0, # where to start the search for a time min
    'nfactor_bm_tail': .9,  # upper bound to start exclude for size
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'], # list of cell types to include
    'quals':[1,2,3], # list of qualities to include
    'plot': True, # boolean for Criticality, if True it will output plots
    'base_saveloc': f'/media/HlabShare/clayton_sahara_work/criticality/', # where to save files. sub directories will be made in this directory
    'verbose':False, # if True criticality will print many statements while running, good for debugging
    'timeout':5000, # in seconds, how long criticality should spend running before moving on to the next block and assuming it won't pass a p-val test
    'none_fact':40, # if bm or tm is None, they're calculated with np.max(burst)/none_fact - the smaller the number the more data EXCLUDE will look through before determining a min and max
    'exclude':True, # boolean, if True then criticality will skip over p-value test if bm or tm surpasses a threshold set by you. meant to save time by not running p-val test if the block isn't going to pass
    'exclude_burst':50, # threshold for excluding size, only matters if exlude is True
    'exclude_time':20, # threshold for excluding size
    'exclude_diff_b':20, # difference between burst min and burst max that if crossed the block is excluded and p_val test is skipped (if the difference between bmin and bmax is < 5, we don't want to use that data)
    'exclude_diff_t':10, # same as diff_b but for time
    'fr_cutoff':50, # if a cell has a FR above this number it's not included
    'start': None, # what block do you want to start at. (1 would start at hour 4 if hour_bins is 4)
    'end': None, # what block do you want to end at. 
    'shuffle':True, # do you want shuffled data to be generated at the end and added to the object?
    'subsample':False, # if True, a subsample of FACTOR cells is pulled ITER times and the AVs are compiled before criticality is run
    'subsample_factor':None, # number of cells pulled
    'subsample_iter':None, # number of times FACTOR cells are pulled
    'subsample_replace':False # if True then 1 cell can be pulled in multiple iterations, if False each cell will only be pulled once.
}
```
*Save*: Boolean - if True then the objects will get saved, if False they will only be returned

*Overlap*: was there an overlap in your data between paths? if so, put the number of hours here. ex: overlap=1

*Timeout*: In seconds - how long criticality should spend before moving onto the next block (necessary because sometimes the p-value test will get stuck and stop progressing and a timeout needs to cut it off and assume it won't pass)


**Outputs**:

*all_objects* : A list of all the criticality objects returned for all paths

*errors* : A list of errors that occured during running. Formatted as [(file, erorr)]

---

### How to run this!

```python
import sahara_work as saw

params = {
    'flag': 2, 
    'ava_binsz': 0.04, 
    'hour_bins': 4,
    'perc': 0.35,
    'bm':None,
    'tm':None,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .9,  
    'nfactor_tm_tail': .9,
    'cell_type': ['FS', 'RSU'],
    'quals':[1,2,3],
    'plot': True,
    'base_saveloc': f'/media/HlabShare/clayton_sahara_work/criticality/',
    'verbose':False,
    'timeout':5000,
    'none_fact':40, 
    'exclude':True, 
    'exclude_burst':50,
    'exclude_time':20,
    'exclude_diff_b':20,
    'exclude_diff_t':10,
    'fr_cutoff':50,
    'save':True,
    'start': None,
    'end': None,
    'shuffle':True,
    'subsample':False,
    'subsample_factor':None,
    'subsample_iter':None, 
    'subsample_replace':False
}

paths = ['32_48/caf19_0326_32_48_scrubbed.npy', '48_64/caf19_0326_48_64_scrubbed.npy', '64_80/caf19_0326_64_80_scrubbed.npy', '80_96/caf19_0326_80_96_scrubbed.npy', '96_112/caf19_0326_96_112_scrubbed.npy'] 

all_data = saw.lilo_and_stitch(paths, params, overlap=0, timeout = params['timeout'], save = params['save'])
```

## Crit Class

This class is meant to be able to be shared with other labs and expanded as needed, thus its inputs are very basic and generalizeable (although lengthy)

**__init__**

In order to init a Crit object the only parameter that is completely necessary is *spikewords*, an n-d array of binarized spike times with shape (ncells, nbins). All other parameters have default values and thus it could be run with only this input. All other possible inputs are: (descriptions found under params of lilo_and_stitch)

*perc*: default 0.35

*nfactor_bm* : default 0

*nfactor_tm* = 0

*nfactor_bm_tail*: default 1

*nfactor_tm_tail*: default 1

*bm*: default None

*tm*: default None

*saveloc*: default ''

*pltname*: default ''

*plot*: default True

*none_fact*: default 40

*exclude*: default False

*exclude_burst*: default0

*exclude_time*: default0

*exclude_diff_b*: default 20

*exclude_diff_t*: default 10

*subsample*: default False

*subsample_factor*: default None

*subsample_iter*: default None

*subsample_replace*: default False

**run_crit**

This is the function that actually runs the criticality analysis using the params given in the init. 

Parameters:

*flag*: 1 for DCC and 2 for P_value test and DCC

*verbose*: If True criticality will output lots of stats while running. good for debugging, annoying for the full analysis


**show_plots**

Displays plots if they exist - the 2 CDFs from the p_value test and the scaling relation plot

**plot_raster**

Plots a raster plot of activity from start to end 

Parameters:
*window_start*: default 200, bin to start plot with 

*window_end*: defualt 400, bin to end plot with 

*saveplot*: default False, save or nah 

*show*: default True, show or nah

### How to use
```python
from sahara_work import Crit_class
import musclebeachtools as mbt
import numpy as np

cells = np.load('/path/to/cells', allow_pickle=True)
good_cells = [cell for cell in cells if cell.quality < 4]
spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = 0.04, binarize = 1)

crit = Crit_class(spikewords)
crit.run_crit(flag = 2, verbose=True)
print(crit.dcc)
```

## Crit_hlab

The *exact* same as crit_class, except this is is built for our data and meant to be run in lilo_and_stitch. Therefore it expects whatever function that runs it to set the animal and date and time frame and paths and so on. Check out the init method to see everything it expects to be set. 

Would 100% recomend using this class as the plots it generates are more tailored to our data, but if you don't double check what you set while running criticality (if you're not using lilo_and_stitch) it may error, so be careful. I tried to fail-safe it but you never know.


# Other (important) Helper Functions
-----
**Informative Functions**:
Unless otherwise stated: animal can be either in 'caf22' format or 'CAF00022' format

These are all dependent on you updating these dictionaries with the proper information after doing a surgery
```python
geno = get_genotype(animal = 'caf22')
print(geno)
'te4'

bday = get_birthday(animal = 'caf22')
print(bday)
'2020-02-17 07:30:00'

sex = get_sex(animal = 'caf22')
print(sex)
'F'

regions = get_regions(animal = 'caf22')
print(regions)
['V1', 'CA1']

probe = get_probe(animal = 'caf22', region = 'CA1')
print(probe)
'probe2'

params = get_params(animal = 'caf22')
print(params)
{'nfactor_bm': 5, 
'nfactor_tm': 0, 
'bm': 20, 
'tm': 8, 
'nfactor_bm_tail': 0.8, 
'nfactor_tm_tail': 0.75, 
'quals': [1, 2]}

probe_params = get_params(animal = 'caf22', region= 'V1')
print(probe_params)
{'nfactor_bm': 5, 
'nfactor_tm': 0, 
'bm': 20, 
'tm': 8, 
'nfactor_bm_tail': 0.8, 
'nfactor_tm_tail': 0.75, 
'quals': [1, 2]}

hstype =  get_hstype(animal = 'caf22')
print(hstype)
['EAB50chmap_00','EAB50chmap_00']
```