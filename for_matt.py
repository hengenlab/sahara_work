from sahara_work import Crit
import numpy as np
import musclebeachtools_hlab as mbt


params = {
    'flag': 2, # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04, # in seconds
    'perc': 0.35,
    'nfactor_bm':0, 
    'nfactor_tm':0,
    'nfactor_bm_tail':.8, # upper bound to start exclude for burst
    'nfactor_tm_tail': .8, # upper bound to start exclude for time 
    'plot' : True
    }

path = '' # path to neuron objects
cells = np.load(path, allow_pickle=True)
good_cells = [cell for cell in cells if cell.quality in [1,3]]
pltname='testing'
saveloc='/media/HlabShare/clayton_sahara_work/criticality/'


spikewords = mbt.spiketimes_to_spikewords(good_cells, binsz=.04, binarize=1)
crit = Crit(spikewords, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                                nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = saveloc,
                                pltname=pltname, plot = params['plot'])

crit.run_crit(flag = params['flag'])

dcc = crit.dcc
p_value_burst = crit.p_value_burst
p_value_t = crit.p_value_t

