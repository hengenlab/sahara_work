## Criticality Tools
These are functions to make running large amounts of data through criticality easier and more efficient 

*SKIP DOWN TO THE LAST FUNCTION!* these are just helper functions

### break_up_mat
**Parameters**:
FR_mat: The binarized spike matrix
small_bin: the binsize used to binarize the matrix
hour_bins: how you want the data broken up

**Output**:
Returns an array of the binarized data reshaped into chunks. So if the data was 12 hours and you wanted it in 3 hour bins, you'll get an array of shape (4,#cells,#bins)

### looped_crit
**Parameters**:
FR_mat: the array of binarized spikes
params: dictionary of parameters
plot: a boolean, if you want the function to output and save plots. this should be false at the moment, i'm making it better

**Outputs**:
Outputs a dictionary with all the possible criticality data you could need


---
## lilo_and_stitch
This bad boy takes in a list of paths, in addition to a dictionary of parameters, and runs criticality for all the paths and saves the master_dictionary to each directory you give it

**Parameters**:


*Paths*: a list of paths to each of the scrubbed neurons you want to run criticality on. Make sure the path is **complete**, and goes all the way to a .npy file. 
ex: paths = ['/media/Hlabshare/clayton_sahara_work/clustering/caf19/0326/0_12/caf_scrubbed.npy', ...]


*Params*:
A dictionary with the necessary parameters, you can put really whatevery you want in here. If you want it saved in the final output, you can throw it in here. This dictionary gets saved to the final output so we can see what our parameters were looking back. But, there are a few necessary inputs so the function can run. The params dictionary needs look *at least* like this:

params = {
    'ava_binsz': 0.04,
    'hour_bins': 4,
    'perc': 0.30,
    'burstM': 15,
    'tM': 5,
    'quality': [1],
    'time_frame': '0326_32-48',
    'animal': 'caf19',
    'notes': ''
}

Before you run the function, define this with the correct parameters. The 'time_frame' element is defined in the loop,so dont worry about that part being correct for each path. 


Overlap: was there an overlap in your data between paths? if so, put the number of hours here. ex: overlap=1


plot: False. not done yet

**Outputs**:
This outputs an array of all the dictionaries. its huge. don't save it. but this way you can look at it after the function runs. 

---

### How to run this!

```python
import Criticality as cr 
import musclebeachtools_hlab.musclebeachtools as mbt 
import numpy as np 
import matplotlib.pyplot as plt 

params = {
    'ava_binsz': 0.04,
    'hour_bins': 4,
    'perc': 0.30,
    'burstM': 15,
    'tM': 5,
    'quality': [1],
    'time_frame': '0326_32-48',
    'animal': 'caf19',
    'notes': ''
}

paths = ['32_48/caf19_0326_32_48_scrubbed.npy', '48_64/caf19_0326_48_64_scrubbed.npy', '64_80/caf19_0326_64_80_scrubbed.npy', '80_96/caf19_0326_80_96_scrubbed.npy', '96_112/caf19_0326_96_112_scrubbed.npy'] 

all_data = lilo_and_stitch(paths, params, overlap=0, plot=False)
```

As this code runs a lot of things will be outputed, at the end of each block finishing itll print out both p_values and the dcc values so you can see it. 


"master_dict" will get saved to each directory. Inside each of these dictionaries has all the criticality data from each step of the process. 

some important elements of master_dict:
all_p_values_t - this will give you all the duration p_values for each block

all_p_values_burst - all the size p_values for each block

all_dcc_values - the dcc values for each block

dcc_ax_block* - the figure returned by running the flag=3 (dcc test) in the criticality code. this is the PDF plot.




