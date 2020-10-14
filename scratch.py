from sahara_work import Crit
import numpy as np 
import pandas as pd
import signal

def signal_handler(signum, frame):
    print("timeout")
    raise Exception('timeout')


caf22_files = ['/media/HlabShare/clayton_sahara_work/criticality/caf22/0512/probe2/Crit_caf22_probe2_0512_0_8_4hrs_perc35_binsz40ms_q1_3_cellsFS_RSU_0.npy', '/media/HlabShare/clayton_sahara_work/criticality/caf22/0512/probe2/Crit_caf22_probe2_0512_0_8_4hrs_perc35_binsz40ms_q1_3_cellsFS_RSU_1.npy']                                                                                                                                                        

caf40_files = ['/media/HlabShare/clayton_sahara_work/criticality/caf40/0914/Crit_caf40_0914_0_12_4hrs_perc35_binsz40ms_q1_2_3_cellsFS_RSU_0.npy','/media/HlabShare/clayton_sahara_work/criticality/caf40/0914/Crit_caf40_0914_0_12_4hrs_perc35_binsz40ms_q1_2_3_cellsFS_RSU_1.npy']          
caf22_objs = []
caf40_objs = []

for f in caf22_files:
    o = np.load(f, allow_pickle=True)[0]
    caf22_objs.append(o)

for f in caf40_files:
    o = np.load(f, allow_pickle=True)[0]
    caf40_objs.append(o)

errors = []
caf22_results = [[],[]]
caf40_results = [[],[]]

all_bins = np.arange(30,510, 10)
all_perc = np.arange(25, 50, 5)
for i, crit in enumerate(caf22_objs):
    crit.saveloc = f'/media/HlabShare/clayton_sahara_work/criticality/param_testing/caf22_{i}'
    for ip, p in enumerate(all_perc):
        print(f"\n\n--SETTING PERC -- {p}\n\n")
        crit.perc = p
        for ib, b in enumerate(all_bins):
            print(f"\n\n--SETTING BINSIZE -- {b}\n\n")
            crit.ava_binsize = b 

            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(600)
            noerr = True

            try:
                crit.run_crit_from_start(save=False)
            except Exception:
                    print('TIMEOUT or ERROR', flush = True)
                    errors.append(f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- ERRORED')
                    noerr=False
            signal.alarm(0)
            if noerr:
                caf22_results[i].append([crit.animal, crit.p_value_burst, crit.p_value_t, crit.dcc, p, b])
            else:
                caf22_results[i].append([crit.animal, np.nan, np.nan, np.nan, p, b])

np.save('param_testing/all_params_caf22_8hours.npy', caf22_results)


for i, crit in enumerate(caf40_objs):
    crit.saveloc = f'/media/HlabShare/clayton_sahara_work/criticality/param_testing/caf40_{i}'
    for ip, p in enumerate(all_perc):
        print(f"\n\n--SETTING PERC -- {p}\n\n")
        crit.perc = p
        for ib, b in enumerate(all_bins):
            print(f"\n\n--SETTING BINSIZE -- {b}\n\n", flush=True)
            crit.ava_binsize = b 

            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(600)
            noerr = True

            try:
                crit.run_crit_from_start(save=False)
            except Exception:
                    print('TIMEOUT or ERROR', flush = True)
                    errors.append(f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- ERRORED')
                    noerr=False
            signal.alarm(0)
            if noerr:
                caf40_results[i].append([crit.animal, crit.p_value_burst, crit.p_value_t, crit.dcc, p, b])
            else:
                caf40_results[i].append([crit.animal, np.nan, np.nan, np.nan, p, b])
np.save('param_testing/all_params_caf40_8hours.npy', caf40_results)