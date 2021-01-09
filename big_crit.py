import numpy as np
import glob
import pandas as pd
from sahara_work import Crit
from sahara_work import lilo_and_stitch
import sahara_work as sw
from datetime import datetime as dt
import signal
import sys
import os


def run(animal = '', probe = '', rerun = True, redo = False):
    s = f'/media/HlabShare/clayton_sahara_work/clustering/{animal}*/*/{probe}*/co/*scored_*.npy'
    print(s)
    og = [f for f in glob.glob(s)]
    print(f'total # of paths: {len(og)}', flush = True)

    csvloc = '/media/HlabShare/clayton_sahara_work/criticality/all_results.csv'

    if redo:
        paths = og
        sw.write_csv_header(csvloc)
        loaded = []
        np.save('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy', loaded)
    else:
        loaded = np.load('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy')
        paths = []
        for p in og:
            if p not in loaded:
                paths.append(p)
    paths = sorted(paths)
    now = dt.now()
    with open('/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt', 'a+') as f:
        f.write('------ JOB START -------- ')
        f.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
        f.write(f'{len(paths)} PATHS TO DO - of this job - lol fuck u right?')

    print(f'Number of paths left to analyze: {len(paths)}', flush = True)

    params = {
        'redo_paths': True, # eventually this should just be deleted entirely but for now, here we are
        'flag': 2,  # 1 is DCC 2 is p_val and DCC
        'ava_binsz': 0.04,  # in seconds
        'hour_bins': 4,  # durration of block to look at
        'perc': 0.35,
        'nfactor_bm': 0,
        'nfactor_tm': 0,
        'nfactor_bm_tail': .9,  # upper bound to start exclude for burst
        'nfactor_tm_tail': .9,  # upper bound to start exclude for time
        'cell_type': ['FS', 'RSU'],
        'plot': True
    }

    bins = np.arange(0, len(paths), 10)
    for i, b in enumerate(bins):
        print(f"\n\n{b} ---- PATHS COMPLETE \n\n", flush = True)
        if b == bins[-1]:
            ps = paths[b:]
        else:
            ps = paths[b: bins[i + 1]]

        all_objs, errors = lilo_and_stitch(ps, params, rerun = rerun)

        results = []
        for o in all_objs:
            results.append([o.animal, o.scored_by, o.probe, o.date, o.time_frame, o.block_num, o.p_value_burst, o.p_value_t, o.dcc, (o.p_value_burst > 0.05 and o.p_value_t > 0.05)])
            err, appended = sw.write_to_results_csv(o, csvloc)

            if err:
                print('something weird happened, this should not have errored')
            else:
                new_path = o.filename
                loaded = np.load('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy')
                loaded = np.append(loaded, new_path)
                np.save('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy', loaded)

        df = pd.DataFrame(results, columns = ['animal', 'scored', 'probe', 'date', 'time_frame', 'block_num', 'p_value_burst', 'p_value_t', 'dcc', 'passed'])

        group = df.groupby(['animal', 'probe', 'date', 'scored'])

        strs = []
        for i, row in group:
            num_passed = row[row["passed"]].count()['passed']
            total_num = row.count()['passed']
            avg_dcc = row.mean()['dcc']
            animal = row['animal'].to_numpy()[0]
            date = row['date'].to_numpy()[0]
            probe = row['probe'].to_numpy()[0]
            scored = row['scored'].to_numpy()[0]
            s = f'{str(animal)} -- {probe} -- {date} -- {scored} -- passed {num_passed}/{total_num} -- avg dcc {avg_dcc}'
            strs.append(s)

        now = dt.now()
        if len(all_objs) > 0:
            with open('/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt', 'a+') as f:
                f.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
                f.write(f'{b} PATHS DONE - of this job')
                for s in strs:
                    f.write(f'{s}\n')
                f.write('\tERRORS:\n')
                for e in errors:
                    f.write(f'\t{e}\n')


if __name__ == "__main__":
    l = len(sys.argv)
    if l > 1:
        animal = sys.argv[1]
        probe = sys.argv[2]
        rerun = True if sys.argv[3] == 'True' else False
        redo = True if sys.argv[4] == 'True' else False
        print(f'specifying animal -- {animal}')
        print(f'specifying probe -- {probe}')
        print(f'rerun -- {rerun}')
        print(f'redo paths -- {redo}')
        run(animal, probe, rerun, redo)
    else:
        run()
