from sahara_work.crit_class import Crit

class Crit_hlab(Crit):

    def __init__(self, **kwargs):
        super().__init__(
            spikewords = kwargs.get('spikewords'),
            perc = kwargs.get('perc'),
            nfactor_bm = kwargs.get('nfactor_bm'),
            nfactor_tm = kwargs.get('nfactor_tm'),
            nfactor_bm_tail = kwargs.get('nfactor_bm_tail'),
            nfactor_tm_tail = kwargs.get('nfactor_tm_tail'),
            saveloc = kwargs.get('saveloc'),
            pltname = kwargs.get('pltname'),
            plot = kwargs.get('plot')
        )

        # all parameters set by run_crit
        self.burst = kwargs.get('burst')
        self.T = kwargs.get('T')
        self.bm = kwargs.get('bm')
        self.tm = kwargs.get('tm')
        self.p_value_burst = kwargs.get('p_value_burst')
        self.p_value_t = kwargs.get('p_value_t')
        self.dcc = kwargs.get('dcc')
        self.scaling_plot = kwargs.get('scaling_plot')
        self.burst_cdf_plot = kwargs.get('burst_cdf_plot')
        self.t_cdf_plot = kwargs.get('t_cdf_plot')
        self.pre = kwargs.get('pre')
        self.fit = kwargs.get('fit')
        self.xmin = kwargs.get('xmin')
        self.xmax = kwargs.get('xmax')
        self.tmin = kwargs.get('tmin')
        self.tmax = kwargs.get('tmax')
        self.alpha = kwargs.get('alpha')


        # all parameters that can be set by lilo and stitch
        self.time_frame = kwargs.get('time_frame')
        self.block_num = kwargs.get('block_num')
        self.qualities = kwargs.get('qualities')
        self.cell_types = kwargs.get('cell_types')
        self.hour_bins = kwargs.get('hour_bins')
        self.ava_binsize = kwargs.get('ava_binsize')
        self.animal = kwargs.get('animal')
        self.date = kwargs.get('date')
        self.final = kwargs.get('final')
        self.cells = kwargs.get('cells')
        self.probe = kwargs.get('probe')

    def __repr__(self):
        print(f'CRIT_HLAB: {self.animal} -- {self.date} -- {self.time_frame}')


