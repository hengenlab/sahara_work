import musclebeachtools_hlab as mbt
import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import seaborn as sns
import numpy as np
import sys
from datetime import datetime
import time
import glob, os
import sahara_work as sw
import pickle
import pandas as pd 
import te4badbrains.te4_paper_utils as te4

if sys.platform == 'darwin':
    # Mac OS
    basedir = 'Volumes'
else:
    # Rig 221
    basedir = 'media' 

if sys.platform == "darwin":
    # mpl.use('TkAgg')
    plt.switch_backend('TkAgg')
else:
    # mpl.use('Agg')
    plt.switch_backend('Agg')






 
plt_dir = f'/{basedir}/HlabShare/AD_paper/figures/figure1/'
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Loop through all of the data files (this takes some time), and count the number of single units in each recording block.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
animals = { 
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
            'caf60': 'te4',
            'caf61': 'e4',
            'caf62': 'te4',
            'caf69': 'wt',
            'caf72': 'te4',
            'caf77': 'wt',
            'caf78': 'te4',
            'caf79': 'te4',
            'caf80': 'te4',
            'caf81': 'wt',
            'caf82': 'wt',
            'caf84': 'te4',
            'caf88': 'wt',
            'caf89': 'wt',
            'caf90': 'wt',
            'caf92': 'wt',
            'caf94': 'wt',
            'caf95': 'wt',
            'caf96': 'wt',
            'caf97': 'wt',
            'eab47': 'te4',
            'eab50': 'wt'
            }

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Run this to generate the neuron count data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
bigdf = pd.DataFrame(columns = ['animal', 'time', 'time_days', 'neurons', 'genotype'])

f_list = {}
for a in animals:

    # if os.path.exists(f'{plt_dir}{a}dataframe_cellcounts_by_time.pkl'):
    #     os.remove(f'{plt_dir}{a}dataframe_cellcounts_by_time.pkl')
    #     print(f'Deleted extant dataframe {a}dataframe_cellcounts_by_time.pkl')
    # else:
    #     pass

    animalid = te4.oldnamenewname(a)
    df = pd.DataFrame(columns = ['animal', 'time', 'neurons', 'genotype'])
    #
    print(f'Starting {a}.')
    probe = sw.get_probe(a,'CA1')

    if probe == -1:
        probe = sw.get_probe(a,'CA1_DG')
    else:
        pass

    f_list[a] = sorted(glob.glob(f'/{basedir}/HlabShare/Clustering_Data/{animalid}/*/*/{probe}/co/*neurons_group0.npy')) 

    for c,f in enumerate(f_list[a]):
        
        print(f'File {c} of {len(f_list[a])}.')
        # Load the data
        try:
            dat = np.load(f, allow_pickle=True)
            dat = te4.kill3s(dat)
            tage = sw.get_age_sec(dat[0].rstart_time, sw.get_birthday(a))
            #g = [True for i in dat if i.quality<3]
            tmp = {'animal': a, 'time' : tage, 'neurons' : len(dat), 'genotype' : animals[a]}
            if len(dat)>150:
                print(f'File: {f}')
                print(f'n neurons is {len(dat)}')
            df = df.append(tmp, ignore_index = True)
            del(dat)
        except Exception as err:
            print(f'Problem with file:{f}')
            print(err)

    df['time_days'] = df['time']/(24*3600)
    bigdf = bigdf.append(df)
    # df_fn = f'{a}dataframe_cellcounts_by_time.pkl'
    # df.to_pickle(plt_dir + df_fn) 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                      Load the dataframes for each animal

animaldfs = glob.glob(plt_dir+'*pkl')
df = pd.DataFrame()
for anml in animaldfs:
    df = df.append(pd.read_pickle(anml))
#-------------------------------------------------------------------------------
#               Figure out the ordering of the onset of the recordings

def addanimal(anid, start_time):
    if anid in animalorder:
        print (f'{anid} is already included.')
    else:
        animalorder[anid] = start_time

animalorder = dict()
ymax = np.zeros(np.unique(df['animal']).size)
for count,aa in enumerate(np.unique(df['animal'])):
    addanimal(aa,np.min(df.loc[df['animal'] == aa]['time_days']))
    ymax[count] = np.max(df.loc[df['animal'] == aa]['neurons'])

sorted_animal_order = dict(sorted(animalorder.items(), key=lambda item: item[1]))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                                       Plot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
fig1, ax1 = plt.subplots(ncols = 1, nrows=np.unique(df['animal']).size, figsize = [20, 5], sharex = True)

for count,aa in enumerate(sorted_animal_order):

    gt = np.unique(df.loc[df['animal'] == aa]['genotype'])[0]
    if gt == 'te4':
        clr = '#0895cb' # purple/pink
    elif gt == 'wt':
        clr = '#08eecb' # purple
    elif gt == 'e4':
        clr = 'xkcd:slate grey' # dark purple


    ax1[count].bar(df.loc[df['animal'] == aa]['time_days'], df.loc[df['animal'] == aa]['neurons'], color = clr)

    #ax1[count].set_ylim([0,100])
    ax1[count].set_ylim([0,ymax.max()])

    if count == np.unique(df['animal']).size-1:
        for axis in ['top','right']:
            ax1[count].spines[axis].set_linewidth(0)
        for axis in ['bottom','left']:
            ax1[count].spines[axis].set_linewidth(1)
        ax1[count].set_xlabel('Time (days)')
    else:
        for axis in ['top','bottom','left','right']:
            ax1[count].spines[axis].set_linewidth(0)
            #ax1[count].tick_params(labelbottom=False)
            ax1[count].set_yticks([])
            ax1[count].tick_params(axis='x', which='both',
                bottom=False)

    ax1[count].set_ylabel(aa, rotation=0)

fig1.suptitle('Number of recorded single units by animal over time.')
#fig1.tight_layout()


today = datetime.today().strftime('%Y-%m-%d')
fig1.savefig(plt_dir + f'timelineplot_{today}.pdf')
