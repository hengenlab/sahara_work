import glob
import numpy as np 
import os
import datetime as dt 

dirs = glob.glob('/media/HlabShare/Sleep_Scoring/CAF00022*/*/*/*/co/') 
print(len(dirs)) 
crashed = [] 
for d in dirs: 
    neurons = glob.glob(d+'*lfp_group0.npy') 
    if len(neurons) < 1: 
        crashed.append(d) 

crashed = [c for c in crashed if 'files_block' not in c]
print(f'{len(crashed)} clustering directories with no neurons_group0 file')
count = [] 
now = dt.datetime.now() 
running = []
for f in crashed: 
    fls = glob.glob(f+'*') 
    if len(fls) > 0: 
        mtimes = [os.path.getmtime(ff) for ff in fls]
        mtime = sorted(mtimes)[-1] 
        t = dt.datetime.fromtimestamp(mtime) 
        elapsed = now - t 
        if elapsed.days < 1:
            running.append(f)
            print('Potentially Currently Clustering ', elapsed, ' ', f) 
        else:
            count.append(f) 
print(len(count), ' folders with files still in them') 

with open('todel.sh', 'w') as shfile:
    for f in count:
        if animal in f and 'yifan' not in f:
            shfile.write('rm '+f+'*\n')


