import numpy as np
import criticality as cr
import musclebeachtools as mbt

class Crit:
    """
    Class to look at criticality stuff

    tbd on init and what not
    """

    def __init__(self, spikewords, perc = 0.35, nfactor_bm = 0, nfactor_tm = 0, nfactor_bm_tail = 1, nfactor_tm_tail = 1, saveloc = '',pltname='', plot = True):

        #required parameters
        self.perc = perc
        self.spikewords = spikewords
        n,b = np.shape(spikewords)
        self.num_cells = n
        self.num_bins = b
        self.saveloc = saveloc
        self.pltname = pltname
        self.plot = plot
        self.nfactor_bm = nfactor_bm
        self.nfactor_tm = nfactor_tm
        self.nfactor_bm_tail = nfactor_bm_tail
        self.nfactor_tm_tail = nfactor_tm_tail

        #all parameters set by run_crit
        self.burst = None
        self.T = None
        self.bm = None
        self.tm = None
        self.p_value_burst = None
        self.p_value_t = None
        self.dcc = None
        self.scaling_plot = None
        self.burst_cdf_plot = None
        self.t_cdf_plot = None
        self.pre = None
        self.fit = None
        self.xmin = None
        self.xmax = None
        self.tmin = None
        self.tmax = None
        self.alpha = None

        #all parameters that can be set by lilo and stitch
        self.time_frame = None
        self.block_num = None
        self.qualities = None
        self.cell_types = None
        self.hour_bins = None
        self.ava_binsize = None
        self.animal = None


    def get_params(self):
        print(f'PARAMETERS:\n'
              f'binsize: {self.ava_binsize}'
              f'perc: {self.perc}'
              f'hour bins: {self.hour_bins}'
              f'qualities: {self.qualities}'
              f'cell types: {self.cell_types}'
              f'nfactor_bm: {self.nfactor_bm}'
              f'nfactor_tm: {self.nfactor_tm}'
              f'nfactor_bm_tail: {self.nfactor_bm_tail}'
              f'nfactor_tm_tail: {self.nfactor_tm_tail}')

    def run_crit(self, flag = 2):
        """
        So this will do the whole calculation. I think i'm gonna break it up into methods even further for re-running purposes

        need to:
        - run get avalanches
        - run av_analyses
        """

        R = cr.get_avalanches(self.spikewords, perc = self.perc)
        burst = R['S']
        T = R['T']
        self.burst = burst
        self.T = T

        self.bm = int(np.max(burst)/20)
        self.tm = int(np.max(T)/20)


        crit_params = {
            'flag' : flag,
            'pltname': self.pltname,
            'saveloc' : self.saveloc,
            'plot' : self.plot,
            'bm' : self.bm,
            'tm' : self.tm
        }

        Result = cr.AV_analysis(burst, T, crit_params, nfactor_bm = self.nfactor_bm, nfactor_tm = self.nfactor_tm, nfactor_bm_tail = self.nfactor_bm_tail,
                             nfactor_tm_tail = self.nfactor_tm_tail)

        if flag == 2:
            self.p_value_burst = Result['P_burst']
            self.p_value_t = Result['P_t']

        self.dcc = Result['df']

        if self.plot:
            self.scaling_plot = Result['scaling_relation_plot']
            if flag == 2:
                self.burst_cdf_plot = Result['burst_cdf']
                self.t_cdf_plot = Result['time_cdf']

        self.pre = Result['pre']
        self.fit = Result['fit']
        self.xmin = Result['xmin']
        self.xmax = Result['xmax']
        self.tmin = Result['tmin']
        self.tmax = Result['tmax']
        self.alpha = Result['alpha']


def __get_totaltime(time_frame):
    start_time = int(time_frame[0:time_frame.find('_')])
    stop_time = int(time_frame[time_frame.find('_')+1:])
    total_time = stop_time-start_time
    return total_time

def __get_paramstr(animal,date, time_frame, hour_bins, perc, ava_binsize, quals, cells, idx):
    qual_str = '_'.join(map(str, quals))
    cell_str = '_'.join(cells)
    s = f'{animal}_{date}_{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsize*1000))}ms_q{qual_str}_cells{cell_str}_{idx}'
    return s

def generate_timeframes(start, end, blocksize):
    ts = np.arange(start, end+1, blocksize)
    time_frames = [str(a)+"_"+str(b) for a, b in zip(ts[:-1], ts[1:])]
    return time_frames

params = {
    'flag': 2, # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04, # in seconds
    'hour_bins': 4,# durration of block to look at
    'perc': 0.35,
    'nfactor_bm':0, 
    'nfactor_tm':0,
    'nfactor_bm_tail':.8, # upper bound to start exclude for burst
    'nfactor_tm_tail': .8, # upper bound to start exclude for time
    'quality': [1,3], # all qualities would be [1,2,3]
    'cell_type': ['FS', 'RSU'], 
    'animal' : 'caf22',
    'saveloc' : "/media/HlabShare/clayton_sahara_work/criticality/caf22/0508/",
    'notes': 'using the new class',
    'time_frame_list':['8_16'], 
    'date': "0508",
    'plot' : True
    }

def lilo_and_stitch(paths, params):
    all_objs = []
    for idx, path in enumerate(paths):
        print(f'Working on ---- {path}')
        time_frame = params['time_frame_list'][idx]
        total_time = __get_totaltime(time_frame)

        num_bins = int(total_time/params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])

        cells = np.load(path, allow_pickle = True)
        good_cells = [cell for cell in cells if cell.quality in params['quality'] and cell.cell_type in params['cell_type']]

        spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1)

        for idx in np.arange(0, num_bins):
            if idx == num_bins - 1:
                data = spikewords[:, (idx * bin_len):]
            else:
                data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]

            param_str = __get_paramstr(params['animal'],params['date'], time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], params['quality'], params['cell_type'], idx)
            crit = Crit(data, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                        nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = params['saveloc'],
                        pltname=param_str, plot = params['plot'])

            crit.run_crit(flag = params['flag'])
            crit.time_frame = time_frame
            crit.block_num = idx
            crit.qualities = params['quality']
            crit.cell_types = params['cell_type']
            crit.hour_bins = params['hour_bins']
            crit.ava_binsize = params['ava_binsz']
            crit.animal = params['animal']
            crit.date = params['date']

            print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}')
            np.save(f'{params["saveloc"]}Crit_{param_str}', crit)
            all_objs.append(crit)

    return all_objs

