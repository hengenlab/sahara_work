from sahara_work import Crit
import numpy as np 

path = '/media/HlabShare/clayton_sahara_work/criticality/caf37/0826/Crit_caf37_0826_16_24_4hrs_perc35_binsz40ms_q1_2_3_cellsFS_RSU_0.npy'
crit = np.load(path, allow_pickle = True)[0]

#2 basic features
cells = crit.cells
spikewords = crit.spikewords