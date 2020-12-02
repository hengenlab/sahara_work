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
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num','bday','rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']

    if re_load:
        with open(csvloc, 'w', newline='') as c:
            w =  csv.DictWriter(c, fieldnames=cols)
            w.writeheader()
        loaded = np.array([])
    else:
        loaded = np.load(loaded_file)
        print(f'Number of paths already loaded: {len(loaded)}')

    count = 0
    for i, p in enumerate(paths):
        if i%5 == 0:
            print(f'#paths: {i}', flush = True)

        if count==100:
            print('saving progress', flush = True)
            np.save(loaded_file, loaded)
            count=0

        if p not in loaded or re_load:
            count+=1
            err, dat = s.write_to_results_csv()
            if err:
                print(to_append)
                errs.append(to_append)
            else:
                s.write_to_csv(to_append, cols, csvloc)
                loaded = np.append(loaded, p)

            if 'LOADED' in p:
                base = p[:p.find('_LOADED')]
                os.rename(p, base+'.npy')

    print('saving final progress', flush = True)
    np.save(loaded_file, loaded)
            
    return errs, loaded