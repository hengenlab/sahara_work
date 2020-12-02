import glob
import sahara_work as s
import csv
import os
import numpy as np


def get_all_results(save, csvloc, loaded_file, re_load):
    '''
    because python refuses to play nice with memory. This will load all crit info into 
    a csv to be loaded into a pandas df at a later time

    it doesn't really work. Gonna eventually re-write this correctly but for now - a bandaid
    '''
    paths = glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/*/*/*/Crit*')
    print(f'Total # of paths: {len(paths)}')
    errs = []
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    if re_load:
        with open(csvloc, 'w', newline='') as c:
            w =  csv.DictWriter(c, fieldnames=cols)
            w.writeheader()
        loaded = []
    else:
        loaded = np.load(loaded_file)
        print(f'Number of paths already loaded: {len(loaded)}')

    for i, p in enumerate(paths):
        if i%5 == 0:
            print(f'#paths: {i}', flush = True)

        if i%100 == 0 and i!=0:
            print('saving progress', flush = True)
            np.save(loaded_file, loaded)

        if p not in loaded or re_load:
            err, to_append = s.lil_helper_boi(p)
            if err:
                print(to_append)
                errs.append(to_append)
            else:
                s.write_to_csv(to_append, cols, csvloc)
                loaded.append(p)

            if 'LOADED' in p:
                base = p[:p.find('_LOADED')]
                os.rename(p, base+'.npy')
            
    return errs, loaded