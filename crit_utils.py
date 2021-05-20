import glob
import sahara_work as s
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns
import criticality as cr
import musclebeachtools as mbt
import csv
import os
import numpy as np
import sahara_work as saw
from datetime import datetime as dt
from datetime import timedelta
from sahara_work import Crit
from sahara_work.crit_hlab import Crit_hlab
import re
import pandas as pd
import os
import signal
import gc
from copy import deepcopy as cdc
import time

def get_cols():
    '''
    Returns the columns for the large data frame I use to look at data
    nice to have them in one place so I don't have to change it in 4 different
    functions when i add a column
    '''
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'scored', 'bday', 'rstart_time', 'age', 'geno',
             'p_val_b', 'p_val_t', 'dcc', 'alpha', 'beta', 'fit_sigma', 'sigma', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t',
              'xmin', 'xmax', 'tmin', 'tmax', 'burstperc', 'Tperc', 'k3b', 'k3t', 'excluded_b', 'excluded_t', 'acc1', 'acc2', 'br1', 'br2', 'burst', 'T', 'num_cells', 'npy_file']   
    return cols

def get_all_results(csvloc, loaded_file, re_load):
    """
    because python refuses to play nice with memory. This will load all crit info into
    a csv to be loaded into a pandas df at a later time

    it doesn't really work. Gonna eventually re-write this correctly but for now - a bandaid


    **** deprecated **** don't use this unless you really know what you're doing and want to
    """
    paths = sorted(glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/*/*/*/Crit*'))
    print(f'Total # of paths: {len(paths)}')
    error_save = '/media/HlabShare/clayton_sahara_work/criticality/results_errors_TODEL.npy'
    if re_load:
        write_csv_header(csvloc)
        loaded = np.array([])
        errs = np.array([])
    else:
        loaded = np.load(loaded_file)
        errs = np.load(error_save, allow_pickle = True)
        print(f'Number of paths already loaded: {len(loaded)}')

    count = 0
    for i, p in enumerate(paths):
        if i % 5 == 0:
            print(f'#paths: {i}', flush = True)
            print(f'count: {count}', flush = True)

        if count == 100:
            print('saving progress', flush = True)
            np.save(loaded_file, loaded)
            np.save(error_save, errs)
            count = 0

        if (p not in loaded and p not in errs) or re_load:
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
    '''
    writes information to a csv file

    data: what to write, should be an array
    cols: the columns corresponding to the array
    loc: what csv you're writing to
    '''
    d = dict(zip(cols, data))
    with open(loc, 'a', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = cols)
        w.writerow(d)

def write_to_pkl(data, cols, loc):
    '''
    appends to the pickle file if it exists, otherwise
    it makes a new one
    '''
    if os.path.exists(loc):
        df = pd.read_pickle(loc)
        temp = pd.DataFrame([data], columns = cols)
        df = df.append(temp)
        df.to_pickle(loc, protocol=4)
    else:
        df = pd.DataFrame([data], columns = cols)
        df.to_pickle(loc, protocol=4)

def write_to_results_csv(crit, loc):
    '''
    pulls info from a crit object and writes it to a csv

    crit: crit object
    loc: csv to write to

    returns: errors if they occur, array of data that was written to the file
    '''
    cols = get_cols()
    err, data = s.lil_helper_boi(crit)
    if err:
        print('this path failed, plz just fucking delete it and re-do this path ffs')
        return err, data
    write_to_csv(data, cols, loc)
    return err, data

def write_to_results_pkl(crit, loc):
    cols = get_cols()
    err, data = s.lil_helper_boi(crit)
    if err:
        print('This path failed, this really shouldnt happen idk what to tell ya')
        print(err)
        return err, data
    write_to_pkl(data, cols, loc)
    return err, data


def write_csv_header(loc):
    '''
    writes the columns into a csv so pandas is nicer later on
    will overwrite the csv if it already exists so be ~careful~

    loc: csv location
    '''
    cols = get_cols()

    with open(loc, 'w', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = cols)
        w.writeheader()


def write_to_results_csv_from_path(p, loc):
    '''
    loads a crit object from p and writes data to loc

    p: path to object
    loc: location of csv

    returns: errors if they occur, data written to file
    '''
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


def get_data_perc(burst, xmin, xmax):
    '''
    given a distribution of AVs (size or durration) and a min and max
    determins the percentage of data that was included in analysis

    returns: percentage (0-1) of data between min and max
    '''
    burst = np.array(burst)
    good_index = np.where(np.logical_and(burst>xmin, burst<xmax))[0]
    perc = len(good_index)/len(burst)
    return perc

def lil_helper_boi(crit):
    '''
    my lil helper boi function

    pulls all necessary info from a crit object to be written to the results csv
    '''
    err = False

    try:
        birth = s.get_birthday(crit.animal)
        start_time = crit.rstart_time
        start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
        age = start_time - birth
        age = age + timedelta(hours = int((crit.block_num * crit.hour_bins)))
        geno = s.get_genotype(crit.animal)
        burstperc = get_data_perc(crit.burst, crit.xmin, crit.xmax)
        Tperc = get_data_perc(crit.T, crit.tmin, crit.tmax)
        if crit.p_value_burst is None or crit.p_value_t is None:
            passed = None
        else:
            passed = (crit.p_value_burst > 0.05 and crit.p_value_t > 0.05)

        info = [crit.animal, crit.probe, crit.date, crit.time_frame, crit.block_num, crit.scored_by, birth, start_time, age, geno,
                crit.p_value_burst, crit.p_value_t, crit.dcc, crit.alpha, crit.beta, crit.fit_sigma, crit.sigma, passed, crit.kappa_burst, crit.kappa_t, crit.k2b, crit.k2t,
                crit.kprob_b, crit.kprob_t, crit.xmin, crit.xmax, crit.tmin, crit.tmax, burstperc, Tperc, crit.k3b, crit.k3t, crit.EXCLUDED_b, crit.EXCLUDED_t, 
                crit.acc1, crit.acc2, crit.br1, crit.br2, crit.burst, crit.T, crit.num_cells, crit.pathname]
    except Exception as e:
        print(f'error: {e}')
        err = True
        info = [e]
    return err, info

def get_paths(scorer = '', geno=None, animal = '', probe = ''):
    '''
    **** deprecated as fuck ****
    meant to pull all paths that need to be analyzed for criticality
    now things are on the chpc so get w it
    '''
    s = f'/media/HlabShare/Clustering_Data/{animal}*/*/*/{probe}*/co/'
    print(s)
    basepaths = [f for f in glob.glob(s) if 'files_block' not in f]
    print(f'total # of folders: {len(basepaths)}', flush = True)
    og = []
    errors = []
    for f in basepaths:
        neuron_files = glob.glob(f+'*neurons_group0*npy')
        scored_files = glob.glob(f+f'*scored_{scorer}*.npy')
        xgb_files = glob.glob(f+f'*neurons_group0.npy')
        animal, _, _, _ = get_info_from_path(f)
        thisgeno = get_genotype(animal)
        if geno is None or thisgeno in geno:
            if scorer == 'xgb':
                og = np.concatenate([og, xgb_files])
            elif len(scored_files) > 0:
                og = np.concatenate([og, scored_files])
            elif len(neuron_files) == 1 and (scorer=='xgb' or scorer == ''):
                og = np.concatenate([og, neuron_files])
            else:
                errors.append(f)
    return og

def get_birthday(animal):
    '''
    animal: string animal name in format 'abc123'
    returns: datetime object of animals birthday at 7:30am on the day
    '''
    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
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
        'caf66': dt(2020, 3, 22, 7, 30),
        'caf69': dt(2020, 1, 5, 7, 30),
        'caf71': dt(2019, 8, 30, 7, 30),
        'caf72': dt(2020, 11, 24, 7, 30),
        'caf73': dt(2020, 1, 5, 7, 30),
        'caf74': dt(2019, 12, 5, 7, 30),
        'caf75': dt(2020, 1, 5, 7, 30),
        'caf77': dt(2020, 1, 5, 7, 30),
        'caf78': dt(2020, 4, 19, 20, 30),
        'caf79': dt(2020, 4, 19, 20, 30),
        'caf80': dt(2020, 4, 19, 20, 30),
        'caf81': dt(2019, 12, 5, 7, 30),
        'caf82': dt(2019, 12, 5, 7,30),
        'caf84':dt(2020, 8, 10, 7, 30),
        'caf88':dt(2019, 12, 5, 7, 30),
        'caf89':dt(2020, 7, 23, 7, 30),
        'caf90':dt(2021, 2, 24, 7, 30),
        'caf91':dt(2021, 2, 24, 7, 30),
        'caf92':dt(2021, 2, 24, 7, 30),
        'caf94':dt(2021, 1, 30, 7, 30),
        'caf95':dt(2021, 1, 30, 7, 30),
        'caf96':dt(2021, 1, 30, 7, 30),
        'caf97':dt(2021, 1, 30, 7, 30),
        'caf99':dt(2020, 10, 1, 7, 30),
        'caf100':dt(2021, 3, 1, 7, 30),
        'caf101':dt(2021, 3, 1, 7, 30),
        'caf102':dt(2021, 1, 30, 7, 30),
        'caf103':dt(2020, 11, 24, 7, 30),
        'eab52': dt(2019, 4, 19, 7, 30),
        'eab47': dt(2019, 2, 17, 7, 30),
        'eab': dt(2019, 2, 17, 7, 30),
        'eab50': dt(2019, 2, 15, 7, 30),
        'eab40': dt(2018, 12, 5, 7, 30)
    }

    return bdays[animal]

def encode_animal(animal):
    '''
    meant to help lower the size of arrays carrying information, turns animal into 
    number (cause strings are large compared to ints)

    needs to be updated
    '''
    animals = ['caf01', 'caf19', 'caf22', 'caf26', 'caf34', 'caf37', 'caf40', 'caf42', 
                'caf48', 'caf49', 'caf50', 'caf52', 'caf54', 'caf55', 'caf58', 'caf60', 
                'caf61', 'caf62', 'caf66', 'caf69', 'caf71', 'caf72', 'caf73', 'caf74', 
                'caf75','caf77', 'caf78', 'caf79', 'caf80', 'caf81', 'caf82', 'caf84',
                'caf88', 'caf89', 'caf90', 'caf91', 'caf92', 'caf94', 'caf95', 'caf96', 
                'caf97','eab52', 'eab47', 'eab', 'eab50', 'eab40']
    nums = np.arange(len(animals))

    keys = dict(zip(animals, nums))

    return keys[animal]

def decode_animal(num):
    '''
    decodes the number from encode_animal() back into a string
    '''

    animals = ['caf01', 'caf19', 'caf22', 'caf26', 'caf34', 'caf37', 'caf40', 'caf42', 
                'caf48', 'caf49', 'caf50', 'caf52', 'caf54', 'caf55', 'caf58', 'caf60', 
                'caf61', 'caf62', 'caf66', 'caf69', 'caf71', 'caf72', 'caf73', 'caf74',
                 'caf75', 'eab52', 'eab47', 'eab', 'eab50', 'eab40']

    return animals[num]

def get_regions(animal):
    '''
    returns a list of regions that was recorded from that animal
    '''
    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
    regions = {
        'caf01': ['CA1'],
        'caf19': ['CA1'],
        'caf22': ['V1','CA1'],
        'caf26': ['M1','CA1','S1'],
        'caf34': ['S1','M1','CA1','NAc'],
        'caf37': ['CA1'],
        'caf40': ['CA1'],
        'caf42': ['M1','ACad','CA1','RSPv','V1'],
        'caf48': ['CA1'],
        'caf49': ['CA1'],
        'caf50': ['CA1'],
        'caf52': ['CA1'],
        'caf54': ['V1'],
        'caf55': ['V1'],
        'caf58': ['CA1'],
        'caf60': ['CA1'],
        'caf61': ['CA1'],
        'caf62': ['CA1'],
        'caf66': ['RSPv','V1','C_Pu','LGN','Sup_col','S1','CA1','M1'],
        'caf69': ['ACad','RSPv','V1','CA1'],
        'caf71': ['V1','RSPv','CA1','ACad'],
        'caf72': ['CA1'],
        'caf73': ['ACad','CA1','RSPv','V1'],
        'caf74': ['ACad','RSPv','CA1','V1'],
        'caf75': ['ACad','CA1','RSPv','V1'],
        'caf77': ['CA1','RSPv','ACad','V1'],
        'caf78': ['CA1'],
        'caf79': ['CA1'],
        'caf80': ['CA1'],
        'caf81': ['ACad','V1', 'CA1', 'RSPv'],
        'caf82': ['CA1','RSPv','V1','ACad'],
        'caf84': ['CA1'],
        'caf88': ['CA1','ACad','RSPv','V1'],
        'caf89': ['CA1'],
        'caf90': ['CA1'],
        'caf91': ['CA1'],
        'caf92': ['CA1'],
        'caf94': ['CA1'],
        'caf95': ['CA1'],
        'caf96': ['CA1'],
        'caf97': ['CA1'],
        'caf99': ['C_Pu','S1','M1','NAc','LGN','V1','RSPv','Sup_col'],
        'caf100': ['CA1'],
        'caf101': ['CA1'],
        'caf102': ['CA1'],
        'caf103': ['CA1'],
        'eab52': ['CA1','V1'],
        'eab47': ['M1_M2','CA1','V2'],
        'eab': ['M1_M2','CA1','V2'],
        'eab50': ['C_Pu','C_Pu','M1_M2','CA1_DG','CA1_DG','S1','Sup_col','V1_V2'],
        'eab40': ['S1','CA1','M1','M2']
    }
    return regions[animal]


def get_sex(animal):
    '''
    returns the sex of an animal recorded from
    '''
    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
    sex = {
        'caf01': 'M',
        'caf19': 'F',
        'caf22': 'F',
        'caf26': 'F',
        'caf34': 'F',
        'caf37': 'M',
        'caf40': 'F',
        'caf42': 'F',
        'caf48': 'M',
        'caf49': 'M',
        'caf50': 'M',
        'caf52': 'F',
        'caf54': 'M',
        'caf55': 'M',
        'caf58': 'F',
        'caf60': 'F',
        'caf61': 'F',
        'caf62': 'F',
        'caf66': 'M',
        'caf69': 'F',
        'caf71': 'M',
        'caf72': 'F',
        'caf73': 'F',
        'caf74': 'M',
        'caf75': 'M',
        'caf77': 'F',
        'caf78': 'F',
        'caf79': 'F',
        'caf80': 'M',
        'caf81': 'F',
        'caf82': 'M',
        'caf84':'F',
        'caf88': 'M',
        'caf89': 'M',
        'caf90': 'M',
        'caf91': 'F',
        'caf92': 'M',
        'caf94': 'M',
        'caf95': 'M',
        'caf96': 'F',
        'caf97': 'F',
        'caf99': 'F',
        'caf100': 'M',
        'caf101': 'F',
        'caf102': 'M',
        'caf103': 'F',
        'eab52': 'F',
        'eab47': 'M',
        'eab': 'M',
        'eab50': 'F',
        'eab40': 'F'
    }

    return sex[animal]

def get_probe(animal, region): 
    '''
    given an animal and a region, it returns the probe that region is on
    returns: probe, as a string ('probe1') to play nice with our paths

    region must be in the correct format
    '''
    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
    probes = np.array(get_regions(animal)) 
    probenum = np.where(probes == region)[0] 
    if len(probenum)==0: 
        print('That region isnt in this animal. Make sure the region formatting is correct.')
        print(' - '.join(['ACAd', 'ACaD', 'CA1', 'CA1_DG', 'C_Pu', 'LGN', 'M1', 'M1_M2',
       'M2', 'NAc', 'RSPv', 'S1', 'Sup_col', 'V1', 'V1_V2', 'V2']))
        return -1 
    probe = f'probe{probenum[0]+1}' 
    return probe 


def get_params(animal, probe):
    # as of 5/05/21 these params produced an effect on all our data
    # if you need to return to these, they're here, don't delete
    keep_these_params = {
        'caf22': {
                    'nfactor_bm': 5,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':8,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf26': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf34': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf37': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.9,
                    'nfactor_tm_tail':0.9,
                    'quals': [1,2]
                    },
        'caf40': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf42': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 10,
                    'tm':8,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf48': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf49': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.78,
                    'quals': [1,2]
                    },
        'caf50': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf52': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':8,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf58': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf60': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.85,
                    'nfactor_tm_tail':0.85,
                    'quals': [1]
                    },
        'caf61': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':20,
                    'nfactor_bm_tail':0.70,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf62': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf66': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf69': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf72': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf77': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':8,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf78': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf79': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf80': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf81': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':12,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf82': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf84': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf88': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf89': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf90': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf92': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.9,
                    'quals': [1,2]
                    },
        'caf95': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf96': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 40,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1]
                    },
        'eab47': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab50': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab40': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    }
    }
    harsh_params = {
        'caf22': {
                    'nfactor_bm': 5,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':8,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf26': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf34': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf37': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf40': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf42': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 10,
                    'tm':8,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf48': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf49': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.78,
                    'quals': [1,2]
                    },
        'caf50': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf52': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':8,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf58': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf60': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 25,
                    'tm':10,
                    'nfactor_bm_tail':0.70,
                    'nfactor_tm_tail':0.8,
                    'quals': [1]
                    },
        'caf61': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':20,
                    'nfactor_bm_tail':0.70,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf62': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf66': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf69': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf72': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf77': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':8,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf78': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf79': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf80': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf81': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':12,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf82': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf84': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf88': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf89': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf90': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf92': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.9,
                    'quals': [1,2]
                    },
        'caf95': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf96': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1]
                    },
        'eab47': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab50': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab40': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    }
    }

    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
    
    base = {
        'nfactor_bm': 0,
        'nfactor_tm':0,
        'bm': 50,
        'tm':20,
        'nfactor_bm_tail':0.8,
        'nfactor_tm_tail':0.8,
        'quals': [1,2]
    }
    # 
    params = {
        'caf22': {
                    'nfactor_bm': 5,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':8,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf26': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf34': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':6,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf37': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.9,
                    'nfactor_tm_tail':0.9,
                    'quals': [1,2]
                    },
        'caf40': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf42': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 10,
                    'tm':9,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf48': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf49': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.78,
                    'quals': [1,2]
                    },
        'caf50': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.85,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf52': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':8,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf58': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf60': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 40,
                    'tm':15,
                    'nfactor_bm_tail':0.9,
                    'nfactor_tm_tail':0.9,
                    'quals': [1]
                    },
        'caf61': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':20,
                    'nfactor_bm_tail':0.65,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf62': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':1,
                    'nfactor_tm_tail':1,
                    'quals': [1,2]
                    },
        'caf66': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf69': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf72': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf77': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':8,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf78': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf79': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf80': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 40,
                    'tm':10,
                    'nfactor_bm_tail':0.9,
                    'nfactor_tm_tail':0.9,
                    'quals': [1,2]
                    },
        'caf81': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':12,
                    'nfactor_bm_tail':0.75,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf82': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf84': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.75,
                    'quals': [1,2]
                    },
        'caf88': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf89': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':10,
                    'nfactor_bm_tail':0.7,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf90': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 20,
                    'tm':15,
                    'nfactor_bm_tail':0.9,
                    'nfactor_tm_tail':0.9,
                    'quals': [1,2]
                    },
        'caf92': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 15,
                    'tm':9,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.85,
                    'quals': [1,2]
                    },
        'caf95': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 35,
                    'tm':10,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'caf96': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1]
                    },
        'caf97': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab47': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 30,
                    'tm':15,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab50': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    },
        'eab40': {
                    'nfactor_bm': 0,
                    'nfactor_tm':0,
                    'bm': 50,
                    'tm':20,
                    'nfactor_bm_tail':0.8,
                    'nfactor_tm_tail':0.8,
                    'quals': [1,2]
                    }
    }

    probe_params = {
        'caf71': {
            'V1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':10,
                'nfactor_bm_tail':0.85,
                'nfactor_tm_tail':0.85,
                'quals': [1,2]
            },
            'RSPv': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 50,
                'tm':20,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            },
            'ACad': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 20,
                'tm':10,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            }
        }, 
        'caf73': {
            'ACad': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 50,
                'tm':20,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            },
            'RSPv': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 50,
                'tm':20,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            },
            'CA1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.9,
                'nfactor_tm_tail':0.9,
                'quals': [1,2]
            }
        },
        'caf74': {
            'ACad': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 40,
                'tm':15,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            },
            'RSPv': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 50,
                'tm':20,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            },
            'CA1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            }
        },
        'caf75': {
            'V1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            },
            'CA1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            }
        },
        'caf69': {
            'ACad': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            },
            'RSPv': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 20,
                'tm':15,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            }
        },
        'caf81': {
            'V1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            }
        },
        'caf82': {
            'ACad': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':20,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            }
        },
        'caf88': {
            'V1': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 50,
                'tm':15,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            },
            'ACad': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 50,
                'tm':20,
                'nfactor_bm_tail':0.75,
                'nfactor_tm_tail':0.75,
                'quals': [1,2]
            },
            'RSPv': {
                'nfactor_bm': 0,
                'nfactor_tm':0,
                'bm': 30,
                'tm':15,
                'nfactor_bm_tail':0.8,
                'nfactor_tm_tail':0.8,
                'quals': [1,2]
            }
        }
    }

    region = get_regions(animal)[int(probe[-1])-1]

    #pref the probe params
    if animal in probe_params.keys():
        if region in probe_params[animal].keys():
            return probe_params[animal][region]

    # if no probe params but normal params return those
    if animal in params.keys():
        return params[animal]
    
    # otherwise base params it is, thank you for visiting
    return base
    
    

        

def get_genotype(animal):
    '''
    returns genotype of animal
    '''
    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
    genos = {
        'caf01': 'e4',
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
        'caf61': 'e4',
        'caf62': 'te4',
        'caf66': 'wt',
        'caf69': 'wt',
        'caf71': 'app_ps1',
        'caf72': 'te4',
        'caf73': 'app_ps1',
        'caf74': 'app_ps1',
        'caf75': 'app_ps1',
        'caf77': 'wt',
        'caf78': 'te4',
        'caf79': 'te4',
        'caf80': 'te4',
        'caf81': 'wt',
        'caf82': 'wt',
        'caf84':'te4',
        'caf88': 'wt',
        'caf89': 'wt',
        'caf90': 'wt',
        'caf91': 'wt',
        'caf92': 'wt',
        'caf94': 'wt',
        'caf95': 'wt',
        'caf96': 'wt',
        'caf97': 'wt',
        'caf99': 'wt',
        'caf100': 'e4',
        'caf101': 'e4',
        'caf102': 'wt',
        'eab52': 'te4',
        'eab47': 'te4',
        'eab': 'te4',
        'eab50': 'wt',
        'eab40': 'wt'
    }

    return genos[animal]

def get_hstype(animal):
    if len(animal) > 6:
        animal = animal[:3].lower() + str(int(animal[3:]))
    regions = {
        'caf01': ['EAB50chmap_00'],
        'caf19': ['EAB50chmap_00'],
        'caf22': ['EAB50chmap_00','EAB50chmap_00'],
        'caf26': ['EAB50chmap_00','APT_PCB','APT_PCB'],
        'caf34': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf37': ['APT_PCB'],
        'caf40': ['APT_PCB'],
        'caf42': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf48': ['hs64'],
        'caf49': ['hs64'],
        'caf50': ['hs64'],
        'caf52': ['hs64'],
        'caf60': ['APT_PCB'],
        'caf61': ['hs64'],
        'caf62': ['hs64'],
        'caf69': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf71': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf72': ['hs64'],
        'caf73': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf74': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf75': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf77': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf78': ['hs64'],
        'caf79': ['hs64'],
        'caf80': ['hs64'],
        'caf81': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf82': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf84':['hs64'],
        'caf88': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf89': ['hs64'],
        'caf90': ['hs64'],
        'caf91': ['hs64'],
        'caf92': ['hs64'],
        'caf94': ['hs64'],
        'caf95': ['hs64'],
        'caf96': ['hs64'],
        'caf97': ['APT_PCB'],
        'caf99': ['APT_PCB','APT_PCB','APT_PCB','APT_PCB','APT_PCB','APT_PCB','APT_PCB','APT_PCB'],
        'caf100': ['hs64'],
        'caf101': ['hs64'],
        'caf102': ['hs64'],
        'caf103': ['hs64'],
        'eab52': ['EAB50champ_00','EAB50champ_00'],
        'eab47': ['EAB50champ_00','EAB50champ_00','EAB50champ_00'],
        'eab50': ['EAB50champ_00','EAB50champ_00','EAB50champ_00','EAB50champ_00','EAB50champ_00','EAB50champ_00','EAB50champ_00','EAB50champ_00'],
    }

def load_crit(path):
    '''
    loads a crit object
    '''
    return np.load(path, allow_pickle = True)[0]

def update_object(old, save_new = False):
    '''
    cute thought. needs to be fixed, probably doesn't need to exist
    '''
    new_obj = Crit_hlab(

        spikewords = old.spikewords, perc = old.perc, nfactor_bm = old.nfactor_bm,
        nfactor_tm = old.nfactor_tm, nfactor_bm_tail = old.nfactor_bm_tail, nfactor_tm_tail = old.nfactor_tm_tail,
        saveloc = old.saveloc, pltname = old.pltname, plot = old.plot, burst = old.burst, T = old.T, tm = old.tm,
        p_value_burst = old.p_value_burst, p_value_t = old.p_value_t, dcc = old.dcc, scaling_plot = old.scaling_plot,
        burst_cdf_plot = old.burst_cdf_plot, t_cdf_plot = old.t_cdf_plot, pre = old.pre, fit = old.fit,
        xmin = old.xmin, xmax = old.xmax, tmin = old.tmin, tmax = old.tmax, alpha = old.alpha, time_frame = old.time_frame,
        block_num = old.block_num, qualities = old.qualities, cell_types = old.cell_types, hour_bins = old.hour_bins,
        ava_binsize = old.ava_binsize, animal = old.animal, date = old.date, final = old.final, cells = old.cells,
        probe = old.probe, filename = old.filename, pathname = old.pathname, scored_by = old.scored_by
        )

    if save_new:
        np.save(old.filename, [new_obj])

    return new_obj

def __get_totaltime(time_frame):
    '''
    given a string in the format '12_24' it tells you how much time is between those
    two hours
    '''
    start_time = int(time_frame[0:time_frame.find('_')])
    stop_time = int(time_frame[time_frame.find('_') + 1:])
    total_time = stop_time - start_time
    return total_time


def __get_paramstr(animal, probe, date, time_frame, hour_bins, perc, ava_binsize, quals, cells, idx):
    '''
    just does some string stuff to make the final path name
    '''
    qual_str = '_'.join(map(str, quals))
    cell_str = '_'.join(cells)
    s = f'{animal}_{probe}_{date}_{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc * 100))}_binsz{str(int(ava_binsize * 1000))}ms_q{qual_str}_cells{cell_str}_{idx}'
    return s


def generate_timeframes(start, end, blocksize):
    '''
    lmfao
    the good old days
    dont need but dont want to delete
    '''
    ts = np.arange(start, end + 1, blocksize)
    time_frames = [str(a) + "_" + str(b) for a, b in zip(ts[:-1], ts[1:])]
    return time_frames


def gen_timeline():
    locs = ['bs001r/rawdata/', 'bs002r', 'bs003r', 'bs004r', 'bs005r', 'bs006r', 'bs007r', 'bs003r/D1/','bs003r/D1_442b/', 'bs004r/D1/', 'bs005r/D1/', 'bs006r/D1/', 'bs007r/D1/', 'bs007r/D1_442b/', 'bs007r/D1_442a/']
    dat = {}
    for loc in locs:
        print(loc)
        restarts = glob.glob(f'/media/{loc}/*/')
        for folder in restarts:
            animal_pattern = '((caf|eab|CAF|EAB)\d{2,})'
            matches = re.findall(animal_pattern, folder)
            if len(matches) > 0:
                animal = matches[0][0]
                animal_clean = animal[:3].lower() + str(int(animal[3:]))
                try:
                    g = saw.get_genotype(animal_clean)
                except KeyError:
                    print(f'dont have records for {animal_clean} ---- skipping')
                    g=None
                    pass
                if g in ['te4', 'wt', 'e4']:
                    if 'D1' in folder:
                        files = sorted(glob.glob(folder+'*.bin'))
                    else:
                        files = sorted(glob.glob(folder+'*/*.bin'))
                    if len(files) > 0:
                        f1 = files[0]
                        f2 = files[-1]
                        d1 = f1[-23:f1.find('.bin')]
                        d2 = f2[-23:f2.find('.bin')]
                        d1 = dt.strptime(d1, '%Y-%m-%d_%H-%M-%S')
                        d2 = dt.strptime(d2, '%Y-%m-%d_%H-%M-%S') 
                        if animal_clean not in dat.keys():
                            print(f'adding {animal_clean}')
                            dat[animal_clean] = {'min':d1, 'max':d2, 'geno':g}
                        else:
                            if d1 < dat[animal_clean]['min']:
                                dat[animal_clean]['min'] = d1
                            if d2 > dat[animal_clean]['max']:
                                dat[animal_clean]['max'] = d2
    for a in dat.keys():
        bday = saw.get_birthday(a)
        min_age = dat[a]['min'] - bday
        max_age = dat[a]['max'] - bday
        dat[a]['min_age'] = int(min_age.total_seconds()/60/60/24)
        dat[a]['max_age'] = int(max_age.total_seconds()/60/60/24)
    return dat

def save_obj(crit):
    '''
    save shit
    gotta get saved as an actual array so numpy doesn't flip
    '''
    to_save = np.array([crit])
    np.save(f'{crit.saveloc}Crit_{crit.pltname}', to_save)


def signal_handler(signum, frame):
    '''
    stuff for signal
    '''
    print("timeout")
    raise Exception('timeout')


def get_info_from_path(path):
    '''
    if you're using a full and completed path from /media/HlabShare/Clustering_Data/
    this function will pull out the animal, restart date, and probe
    and return all that
    '''
    animal_pattern = '((caf|eab|CAF|EAB)\d{2,})'
    matches = re.findall(animal_pattern, path)
    animal = matches[0][0]

    date_pattern = '\d{8}'
    matches = re.findall(date_pattern, path)
    date = matches[0]
    date = date[0:4]+date[-2:]

    time_frame_pattern = '/\d{1,}_\d{1,}'
    matches = re.findall(time_frame_pattern, path)
    time_frame = matches[0][1:]

    probe = path[path.find('probe'):path.find('probe') + 6]

    animal_clean = animal[:3].lower() + str(int(animal[3:]))

    return animal_clean, date, time_frame, probe

def get_cell_stats(cell):
    '''
    nope
    '''
    fr, xbins = cell.plotFR(binsz = 3600, start = False, end = False,
                            lplot = 0)

    isis = []
    bins = np.arange(cell.start_time, cell.end_time, 300)  # 5 min bins

    for i in range(len(bins)):
        if i == len(bins) - 1:
            end_bin = cell.end_time
        else:
            end_bin = bins[i + 1]
        spk_idxs = np.where(np.logical_and(cell.spike_time_sec > bins[i], cell.spike_time_sec < end_bin))
        isis.append(np.diff(cell.spike_time[spk_idxs]))

    means = np.array([np.mean(i) for i in isis])
    stds = np.array([np.std(i) for i in isis])
    cvs = stds / means
    binned_cvs = [cvs[i * 12:(i + 1) * 12] for i in np.arange(int(xbins[-1]))]
    cv = [np.mean(i) for i in binned_cvs]
    if len(fr) != len(cv):
        print('this isnt working, fix it')
    return fr, cv


def construct_fr_df(paths):
    '''
    **** deprecated **** dont use
    but if you need to canabalize some code - go for it
    '''
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
    seconds_in_day = 60 * 60 * 24
    with open('/media/HlabShare/clayton_sahara_work/criticality/cell_stats.csv', 'w', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = ['animal', 'rstart_time', 'age_start', 'days_old', 'hours_old', 'cell_idx', 'quality', 'fr', 'cv', 'wf'])
        w.writeheader()

    done = []
    final = []
    for i, p in enumerate(paths):
        print(i)
        path_info = p[:p.find('_perc')]
        if path_info in done:
            print('already done')
        else:
            done.append(path_info)
            crit = None
            try:
                crit = np.load(p, allow_pickle = True)[0]
                birth = bdays[crit.animal]
                for cell in crit.cells:
                    start_time = crit.cells[0].rstart_time
                    start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
                    age = start_time - birth

                    fr, cv = get_cell_stats(cell)

                    for i in range(len(fr)):
                        age_now = age + timedelta(hours = i)
                        days_old = age_now.total_seconds() / seconds_in_day
                        hours_old = age_now.total_seconds() / 3600

                        with open('/media/HlabShare/clayton_sahara_work/criticality/cell_stats.csv', 'a', newline = '') as c:
                            w = csv.DictWriter(c, fieldnames = ['animal', 'rstart_time', 'age_start', 'days_old', 'hours_old', 'cell_idx', 'quality', 'fr', 'cv', 'wf'])

                            d = {
                                'animal': crit.animal,
                                'rstart_time': start_time,
                                'age_start': age_now,
                                'days_old': days_old,
                                'hours_old': hours_old,
                                'cell_idx': cell.clust_idx,
                                'quality': cell.quality,
                                'fr': fr[i],
                                'cv': cv[i],
                                'wf': cell.waveform
                            }
                            w.writerow(d)
            except Exception:
                print('This object is final or weird. Go back and find these cells individually to add to the csv')
                final.append(p)
    np.save('/media/HlabShare/clayton_sahara_work/fr_csv_done.npy', done)
    np.save('/media/HlabShare/clayton_sahara_work/fr_csv_FINAL.npy', final)


    if obj.final:
        print('This crit object is final, there are no cells saved here. If youd like to rerun this block start from lilo_and_stitch')
        return
    total_time = __get_totaltime(obj.time_frame)
    num_bins = int(total_time / obj.hour_bins)
    bin_len = int((obj.hour_bins * 3600) / obj.ava_binsize)
    good_cells = [cell for cell in obj.cells if cell.quality in obj.qualities and cell.cell_type in obj.cell_types]
    spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = obj.ava_binsize, binarize = 1)
    idx = obj.block_num
    if idx == num_bins - 1:
        data = spikewords[:, (idx * bin_len):]
    else:
        data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]
    obj.spikewords = data
    param_str = __get_paramstr(obj.animal, obj.probe, obj.date, obj.time_frame, obj.hour_bins, obj.perc, obj.ava_binsize, obj.qualities, obj.cell_types, idx)
    obj.pltname = param_str
    obj.run_crit(flag = flag)
    print(f'BLOCK RESULTS: P_vals - {obj.p_value_burst}   {obj.p_value_t} \n DCC: {obj.dcc}')
    if save:
        to_save = np.array([obj])
        np.save(f'{obj.saveloc}Crit_{param_str}', to_save)
    return obj

def get_age_sec(start_time, birthday):
    '''
    given a start time of a recording (found in the neuron objects)
    and a date time birthday of an animal
    this will return the age of the animal in seconds at the start of the 
    recording
    '''
    start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S') 
    age = start_time - birthday
    seconds = age.total_seconds()
    return seconds

# function names are getting a bit ambiguous
# this shuffles spiketimes and returns shuffled burst and time distributions
def cha_cha_slide(cells, ava_binsize, perc, start = False, return_cells = False):
    '''
    shuffles spiketimes and returns shuffled burst and time distributions
    '''
    shuffled_cells = cdc(cells)
    for cell in shuffled_cells:
        cell.shuffle_times()
    
    spikewords = mbt.n_spiketimes_to_spikewords(shuffled_cells, binsz = ava_binsize, binarize = 1, start = start)
    R = cr.get_avalanches(spikewords, perc = perc)
    burstS = R['S']
    TS = R['T']
    if return_cells:
        return burstS, TS, shuffled_cells
    else:
        return burstS, TS


params = {
    'flag': 2,  # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04,  # in seconds
    'hour_bins': 4,  # durration of block to look at
    'perc': 0.35,
    'bm':None,
    'tm':None,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .9,  # upper bound to start exclude for burst
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time 
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
    'subsample_replace':False,
    'branching_ratio': False,
    'br_binsize': 0.004,
    'br_kmax': 500,
    'br_binlen': 5, # in minutes
    'br_numbins': 3 # begining middle end - this isn't implemented right now, it's just in the middle
}


#this bad boy reruns paths if it fails a pval. not the best, wouldn't recomend
def lilo_and_stitch_extended_edition(paths, params, rerun = False, save = True, overlap = False, verbose = True, timeout = 600):
    all_objs = []
    errors = []
    for idx, path in enumerate(paths):
        tic = time.time()
        basepath = path[:path.rfind('/')]
        
        print(f'\n\nWorking on ---- {path}', flush = True)
        animal, date, time_frame, probe = get_info_from_path(path)
        print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
        total_time = __get_totaltime(time_frame)
        saveloc = os.path.join(params['base_saveloc'], animal, date, probe) + '/'
        print(f'saveloc: {saveloc}', flush=True)
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        if path.find('scored') < 0:
            scorer = 'xgb'
        else:
            scorer = path[path.find('scored')+7:path.find('.npy')]

        num_bins = int(total_time / params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])


        quals = [1, 2, 3]
        fr_cutoff = 50
        try:
            cells = np.load(path, allow_pickle = True)
            good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            num_cells = len(good_cells)
    
            if len(good_cells) < 10:
                quals = [1, 2, 3]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

            elif len(good_cells) < 60:
                quals = [1, 2]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            elif len(good_cells) >= 60:
                quals = [1]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

                # if len(good_cells) > 100:
                #     cell_idxs = np.random.choice(len(good_cells), 50, replace=False)
                #     good_cells = good_cells[cell_idxs]
            if overlap :
                start = 3600
            else:
                start = False
            spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
        except Exception as err:
            print("Neuron File Won't Load")
            print(err)
            errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- ALL --- {scorer} --- ERRORED', path])
            continue
        for idx in np.arange(0, num_bins):
            liltic = time.time()
            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(timeout)
            noerr = True
            try:
                print(f'Working on block {idx} --- hours {idx * params["hour_bins"]}-{(idx + 1) * params["hour_bins"]}', flush = True)
                if idx == num_bins - 1:
                    data = spikewords[:, (idx * bin_len):]
                else:
                    data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]

                param_str = __get_paramstr(animal, probe, date, time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], quals, params['cell_type'], idx)
                crit = Crit_hlab(spikewords = data, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = saveloc,
                            pltname = f'{param_str}_{scorer}', plot = params['plot'])

                crit.run_crit(flag = params['flag'], verbose = verbose)
                crit.time_frame = time_frame
                crit.block_num = idx
                crit.qualities = quals
                crit.cell_types = params['cell_type']
                crit.hour_bins = params['hour_bins']
                crit.ava_binsize = params['ava_binsz']
                crit.animal = animal
                crit.date = date
                crit.final = False
                crit.cells = [cell for cell in cells if cell.quality < 4]
                crit.probe = probe
                crit.scored_by = scorer
                crit.pathname = path
                crit.filename = f'{saveloc}Crit_{param_str}_{scorer}'

            except Exception as err:
                print('TIMEOUT or ERROR', flush = True)
                print(err)
                errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED', path])
                noerr = False
                signal.alarm(0)

            if rerun and noerr:
                while crit.p_value_burst < 0.05 or crit.p_value_t < 0.05:
                    signal.signal(signal.SIGALRM, signal_handler)
                    signal.alarm(timeout)
                    print('\nRERUNNING BLOCK', flush = True)
                    if crit.nfactor_tm_tail < 0.75 or crit.nfactor_bm_tail < 0.75:
                        print('DONE RERUNNNING -- BLOCK WILL NOT PASS\n')
                        signal.alarm(0)
                        break
                    if crit.p_value_burst < 0.05:
                        crit.nfactor_bm_tail -= 0.05
                        #crit.bm += 5
                    if crit.p_value_t < 0.05:
                        crit.nfactor_tm_tail -= 0.05
                        #crit.tm += 5
                    try:
                        crit.run_crit(flag = params['flag'], verbose = verbose)

                    except Exception:
                        print('TIMEOUT or ERROR', flush = True)
                        errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED', path])
                        signal.alarm(0)
                        noerr = False
                        break
                    signal.alarm(0)

            if noerr:
                print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}', flush = True)
                if save:
                    to_save = np.array([crit])
                    np.save(crit.filename, to_save)
                all_objs.append(crit)
            liltoc = time.time()
            print(f'Time for 1 block: {(liltoc-liltic)/60} min')
        toc = time.time()
        print(f'TOTAL PATH TIME: {(toc-tic)/60} min')


    return all_objs, errors


def lilo_and_stitch(paths, params, save = True, overlap = False, timeout = False):
    '''
    read the readme online. too much to put here. i'll even link it for you
    https://github.com/hengenlab/sahara_work/blob/master/README.md
    '''
    all_objs = []
    errors = []
    for i, path in enumerate(paths):
        tic = time.time()
        basepath = path[:path.rfind('/')]
        
        print(f'\n\nWorking on ---- {path}', flush = True)
        animal, date, time_frame, probe = get_info_from_path(path)
        print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
        total_time = __get_totaltime(time_frame)
        saveloc = os.path.join(params['base_saveloc'], animal, date, probe) + '/'
        print(f'saveloc: {saveloc}', flush=True)
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        if path.find('scored') < 0:
            scorer = 'xgb'
        else:
            scorer = path[path.find('scored')+7:path.find('.npy')]
        start_bin = params['start']
        num_bins = params['end']

        if start_bin is None:
            start_bin = 0
        if num_bins is None:
            num_bins = int(total_time / params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])
        print('start: ', start_bin, ' end ', num_bins, ' bin_len ', bin_len)

        # set individual params
        indiv = get_params(animal, probe)
        for key in indiv.keys():
            params[key] = indiv[key]

        quals = params['quals']
        fr_cutoff = params['fr_cutoff']
        try:
            cells = np.load(path, allow_pickle = True)
            good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            
            num_cells = len(good_cells)

            if overlap :
                print('There is an overlap so cutting out the first hour')
                start = 3600
            else:
                start = False
            spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
        except Exception as err:
            print("Neuron File Won't Load")
            print(err)
            errors.append([animal, probe, date, time_frame, 'ALL', scorer, path, err, dt.now(), 'cells'])
            continue

        for idx in np.arange(start_bin, num_bins):
            liltic = time.time()
            if timeout is not False:
                signal.signal(signal.SIGALRM, signal_handler)
                signal.alarm(timeout)
            noerr = True
            try:
                print(f'Working on block {idx} --- hours {idx * params["hour_bins"]}-{(idx + 1) * params["hour_bins"]}', flush = True)
                if idx == num_bins - 1:
                    print('last bin')
                    data = spikewords[:, (idx * bin_len):]
                else:
                    data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]
                print('shape data: ', np.shape(data), ' shape total: ', np.shape(spikewords))
                param_str = __get_paramstr(animal, probe, date, time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], quals, params['cell_type'], idx)
                crit = Crit_hlab(spikewords = data, perc = params['perc'], bm = params['bm'], tm = params['tm'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], none_fact = params['none_fact'], saveloc = saveloc,
                            pltname = f'{param_str}_{scorer}', plot = params['plot'], exclude = params['exclude'], exclude_burst = params['exclude_burst'], exclude_time = params['exclude_time'], 
                            exclude_diff_b = params['exclude_diff_b'], exclude_diff_t=params['exclude_diff_t'], subsample = params['subsample'],
                            subsample_factor = params['subsample_factor'], subsample_iter = params['subsample_iter'], subsample_replace = params['subsample_replace'])
                
                crit.run_crit(flag = params['flag'], verbose = params['verbose'])
                crit.time_frame = time_frame
                crit.block_num = idx
                crit.qualities = quals
                crit.cell_types = params['cell_type']
                crit.hour_bins = params['hour_bins']
                crit.ava_binsize = params['ava_binsz']
                crit.animal = animal
                crit.date = date
                crit.final = False
                crit.rstart_time = good_cells[0].rstart_time
                #crit.all_cells = [cell for cell in cells if cell.quality < 4]
                crit.used_cells = [cell.clust_idx for cell in good_cells]
                crit.probe = probe
                crit.scored_by = scorer
                crit.pathname = path
                crit.filename = f'{saveloc}Crit_{param_str}_{scorer}'
                if params['shuffle']:
                    burstS, TS = cha_cha_slide(good_cells, params['ava_binsz'], params['perc'], start = start, return_cells = False)
                    crit.burstS = burstS
                    crit.TS = TS

            except Exception as err:
                print('TIMEOUT or ERROR', flush = True)
                print(err)
                errors.append([animal, probe, date, time_frame, idx, scorer, path, err, dt.now(), 'crit'])
                noerr = False
                if timeout is not False:
                    signal.alarm(0)
            
            if params['branching_ratio']:
                try:
                    start_hour = idx*params['hour_bins']
                    if idx == num_bins -1:
                        end_hour = int(good_cells[0].end_time/60/60)
                    else:
                        end_hour = (idx+1) * params['hour_bins']
                    mid = (end_hour - start_hour)/2
                    binlen = params['br_binlen']/60
                    s = start_hour + mid
                    e = round(s + binlen, 2)
                    crit.run_branching_ratio(cells = good_cells, binsize=params['br_binsize'], start = s, end = e, kmax = params['br_kmax'])
                    print(f'BRANCHING RATIO: {crit.acc1}, {crit.acc2}, {crit.br1}, {crit.br2}')
                except Exception as err:
                    print('ERROR IN BRANCHING RATIO')
                    print(err)
                    errors.append([animal, probe, date, time_frame, idx, scorer, path, err, dt.now(), 'branching ratio'])
            if noerr:
                print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}', flush = True)
                if save:
                    to_save = np.array([crit])
                    np.save(crit.filename, to_save)
                all_objs.append(crit)
            liltoc = time.time()
            print(f'Time for 1 block: {(liltoc-liltic)/60} min')
        toc = time.time()
        print(f'TOTAL PATH TIME: {(toc-tic)/60} min')


    return all_objs, errors