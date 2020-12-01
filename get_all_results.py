import glob
import sahara_work as s
import csv
import os


def get_all_results(save = True, csvloc = '', re_load=False):
    '''
    because python refuses to play nice with memory. This will load all crit info into 
    a csv to be loaded into a pandas df at a later time

    fingers crossed this works
    '''
    paths = glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/{animal}*/*/{probe}*/Crit*')
    print(f'Total # of paths: {len(paths)}')
    errs = []
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    if re_load:
        with open(csvloc, 'w', newline='') as c:
            w =  csv.DictWriter(c, fieldnames=cols)
            w.writeheader()

    for i, p in enumerate(paths):
        if i%5 == 0:
            print(f'#paths: {i}', flush = True)

        if 'LOADED' not in p or re_load:
            err, to_append = s.lil_helper_boi(p)
            if err:
                print(to_append)
                errs.append(to_append)
            else:
                s.write_to_csv(to_append, cols, csvloc)

            if 'LOADED' in p:
                base = p[:p.find('_LOADED')]
            else:
                base = p[:-4]
            n = base+'_LOADED.npy'
            os.rename(p, n)
            
    return errs