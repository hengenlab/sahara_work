import glob
import sahara_work as s
import csv
import os
import numpy as np
from datetime import datetime as dt
from datetime import timedelta


def get_all_results(csvloc, loaded_file, re_load):
    """
    because python refuses to play nice with memory. This will load all crit info into
    a csv to be loaded into a pandas df at a later time

    it doesn't really work. Gonna eventually re-write this correctly but for now - a bandaid
    """
    paths = glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/*/*/*/Crit*')
    print(f'Total # of paths: {len(paths)}')
    error_save = '/media/HlabShare/clayton_sahara_work/criticality/results_errors_TODEL.npy'
    if re_load:
        cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'bday', 'rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
        with open(csvloc, 'w', newline = '') as c:
            w = csv.DictWriter(c, fieldnames = cols)
            w.writeheader()
        loaded = np.array([])
        errs = np.array([])
    else:
        loaded = np.load(loaded_file)
        errs = np.load(error_save)
        print(f'Number of paths already loaded: {len(loaded)}')

    count = 0
    for i, p in enumerate(paths):
        if i % 5 == 0:
            print(f'#paths: {i}', flush = True)

        if count == 100:
            print('saving progress', flush = True)
            np.save(loaded_file, loaded)
            np.save(error_save, errs)
            count = 0

        if p not in loaded or re_load:
            count += 1
            err, dat = s.write_to_results_csv_from_path(p, csvloc)
            if err:
                if len(dat) < 2:
                    dat = np.append(dat, p)
                errs = np.append(errs, dat)
            else:
                loaded = np.append(loaded, p)

            if 'LOADED' in p:
                base = p[:p.find('_LOADED')]
                os.rename(p, base + '.npy')

    print('saving final progress', flush = True)
    np.save(loaded_file, loaded)
    np.save(error_save, errs)

    return errs, loaded


def write_to_csv(data, cols, loc):
    d = dict(zip(cols, data))
    with open(loc, 'a', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = cols)
        w.writerow(d)


def write_to_results_csv(crit, loc):
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'bday', 'rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    err, data = s.lil_helper_boi(crit)
    if err:
        print('this path failed, plz just fucking delete it and re-do this path ffs')
        return err, data
    write_to_csv(data, cols, loc)
    return err, data


def write_to_results_csv_from_path(p, loc):
    err = False
    crit = None
    try:
        crit = np.load(p, allow_pickle = True)
        crit = crit[0]
    except Exception as er:
        print("won't load object", flush = True)
        err = True
        errors = [p, er]
        return err, errors

    err, data = write_to_results_csv(crit, loc)
    return err, data


def lil_helper_boi(crit):
    err = False

    try:
        birth = s.get_birthday(crit.animal)
        start_time = crit.cells[0].rstart_time
        start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
        age = start_time - birth
        age = age + timedelta(hours = int((crit.block_num * crit.hour_bins)))
        geno = s.get_genotype(crit.animal)
        info = [crit.animal, crit.probe, crit.date, crit.time_frame, crit.block_num, birth, start_time, age, geno, crit.p_value_burst, crit.p_value_t, crit.dcc, (crit.p_value_burst > 0.05 and crit.p_value_t > 0.05), crit.kappa_burst, crit.kappa_t, crit.k2b, crit.k2t, crit.kprob_b, crit.kprob_t]
    except Exception as e:
        print(f'error: {e}')
        err = True
        info = [e]
    return err, info


def get_birthday(animal):
    bdays = {
        'caf01': dt(2019, 12, 24, 7, 30),
        'caf19': dt(2020, 1, 19, 7, 30),
        'caf22': dt(2020, 2, 17, 7, 30),
        'caf26': dt(2020, 2, 20, 7, 30),
        'caf34': dt(2020, 3, 18, 7, 30),
        'caf37': dt(2019, 8, 18, 7, 30),
        'caf40': dt(2020, 2, 20, 7, 30),
        'caf42': dt(2020, 2, 20, 7, 30),
        'caf48': dt(2020, 7, 20, 7, 30),
        'caf49': dt(2020, 7, 20, 7, 30),
        'caf50': dt(2020, 7, 20, 7, 30),
        'caf52': dt(2020, 4, 19, 7, 30),
        'caf54': dt(2020, 7, 11, 7, 30),
        'caf55': dt(2020, 7, 11, 7, 30),
        'caf58': dt(2020, 9, 23, 7, 30),
        'caf60': dt(2020, 9, 23, 7, 30),
        'caf61': dt(2019, 12, 11, 7, 30),
        'caf62': dt(2019, 11, 18, 7, 30),
        'eab52': dt(2019, 4, 19, 7, 30),
        'eab47': dt(2019, 2, 17, 7, 30),
        'eab': dt(2019, 2, 17, 7, 30),
        'eab50': dt(2019, 2, 15, 7, 30),
        'eab40': dt(2018, 12, 5, 7, 30)
    }

    return bdays[animal]


def get_genotype(animal):
    genos = {
        'caf01': 'te4',
        'caf19': 'te4',
        'caf22': 'te4',
        'caf26': 'wt',
        'caf34': 'wt',
        'caf37': 'te4',
        'caf40': 'wt',
        'caf42': 'wt',
        'caf48': 'te4',
        'caf49': 'te4',
        'caf50': 'e4',
        'caf52': 'te4',
        'caf54': 'myt1l',
        'caf55': 'myt1l',
        'caf58': 'e4',
        'caf60': 'te4',
        'caf61': '?',
        'caf62': 'te4',
        'eab52': 'te4',
        'eab47': 'te4',
        'eab': 'te4',
        'eab50': 'wt',
        'eab40': 'wt'
    }

    return genos[animal]
