
params = {
    'flag': 2,  # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04,  # in seconds
    'hour_bins': 4,  # durration of block to look at
    'perc': 0.35,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .9,  # upper bound to start exclude for burst
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'],
    'plot': True,
    'quals': None
}


def lilo_and_stitch(paths, params, rerun = False, save = True, overlap = False):
    all_objs = []
    errors = []
    for idx, path in enumerate(paths):
        basepath = path[:path.rfind('/')]
        
        print(f'\n\nWorking on ---- {path}', flush = True)
        animal, date, time_frame, probe = get_info_from_path(path)
        print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
        total_time = __get_totaltime(time_frame)
        saveloc = f'/media/HlabShare/clayton_sahara_work/criticality/param_testing/{animal}/{date}/{probe}/'
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        if path.find('scored') < 0:
            scorer = 'xgb'
        else:
            scorer = path[path.find('scored')+7:path.find('.npy')]

        num_bins = int(total_time / params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])


        quals = params['quals']
        fr_cutoff = 50
        try:
            cells = np.load(path, allow_pickle = True)
            good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            num_cells = len(good_cells)

            print(f'WITH QUALS: {quals}, NCELLS: {num_cells}')
            spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
        except Exception as err:
            print("Neuron File Won't Load")
            print(err)
            pass
        for idx in np.arange(0, num_bins):
            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(600)
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

                crit.run_crit(flag = params['flag'])
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
                    signal.alarm(900)
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
                        crit.run_crit(flag = params['flag'])

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

        with open(f'{basepath}/done.txt', 'w+') as f:
            f.write('done')

    return all_objs, errors
