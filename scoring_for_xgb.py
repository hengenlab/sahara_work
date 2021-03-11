import numpy as np 
import glob
import sys
import musclebeachtools as mbt

existing = []

with open('/media/HlabShare/models/caf/base.txt', 'r') as f:
    for line in f:
        existing.append(line.strip())


fls = sorted(glob.glob('/media/HlabShare/Clustering_Data/*/*/*/*/co/*neurons_group0.npy'))
rand_ints = np.unique(np.random.randint(0,len(fls),size=10))
model = '/media/HlabShare/models/caf/xgb_model_k10nfrpr_prob2_0_corr10_caf_sah_basedev3'

for i, ri in enumerate(rand_ints):
    keepgoing = input(f"\nOn file {i} of {len(rand_ints)}. Wanna keep going or nah?  ")

    while keepgoing.lower() not in ['y', 'yes', 'sure', 'might as well', 'fuck it', 'ugh', 'n', 'no', 'nah', 'fuck you', 'fuck right off', 'i hate you', 'i hate myself', 'plz no']:
        keepgoing = input('\nok i know youre upset but lets keep it together shall we? wanna do another file? ')

    if keepgoing.lower() in ['y', 'yes', 'sure', 'might as well', 'fuck it', 'ugh']:
        if fls[ri] not in existing:
            fl = fls[ri]
            print(fl)
            cells = np.load(fl, allow_pickle=True)
            mbt.autoqual(cells,model)
            for cell in cells:
                cell.quality = 0
                cell.checkqual()
            np.save(fl,cells)
            with open('/media/HlabShare/models/caf/base.txt', 'a+') as f:
                f.write(fl + '\n')
        else:
            print('\n-----already done this file, yeet-----\n')
            pass
    elif keepgoing.lower() in ['n', 'no', 'nah', 'fuck you', 'fuck right off', 'i hate you', 'i hate myself', 'plz no']:
        print("\nYea I feel that. bye")
        sys.exit()
    else:
        print('\nk i dont know how to take that. so were just gonna do another file\n')


# rescoring 

allfiles = glob.glob('/media/HlabShare/models/model_acc_test/*')
changes = np.load('/media/HlabShare/clayton_sahara_work/criticality/testing/to_change2.npy', allow_pickle=True)
changes = dict(changes)

for i, f in enumerate(allfiles):
    print(f'File: {i} of {len(allfiles)} --- {f}')
    cells = np.load(f, allow_pickle=True)
    name = f[f.find('test')+5:]
    try:
        tofix = changes[f]
    except KeyError:
        tofix = []
    for cell in cells:
        if cell.mean_amplitude < 35:
            print('Low amp - fixing')
            cell.quality = 4
        elif cell.clust_idx in tofix:
            cell.checkqual()
    np.save('/media/HlabShare/models/model_acc_test_updated/' + name, cells)


# fixing cells
fls = sorted(glob.glob('/media/HlabShare/Clustering_Data/*/*/*/*/co/*neurons_group0.npy'))
rand_ints = np.unique(np.random.randint(0,len(fls),size=10))
model = '/media/HlabShare/models/caf/xgb_model_k10nfrpr_prob2_0_corr10_caf_sah_basedev4'

for i, ri in enumerate(rand_ints):
    keepgoing = input(f"\nOn file {i} of {len(rand_ints)}. Wanna keep going or nah?  ")

    while keepgoing.lower() not in ['y', 'yes', 'sure', 'might as well', 'fuck it', 'ugh', 'n', 'no', 'nah', 'fuck you', 'fuck right off', 'i hate you', 'i hate myself', 'plz no']:
        keepgoing = input('\nok i know youre upset but lets keep it together shall we? wanna do another file? ')

    if keepgoing.lower() in ['y', 'yes', 'sure', 'might as well', 'fuck it', 'ugh']:
        tofix = []
        fl = fls[ri]
        print(fl)
        cells = np.load(fl, allow_pickle=True)
        mbt.autoqual(cells,model)
        for cell in cells:
            OGq = cell.quality
            cell.checkqual()
            newQ = cell.quality
            if OGq != newQ:
                print('\nNEW QUAL DIFFERENT - saving cell')
                tofix.append(cell)
        if len(tofix) > 0:
            print('---saving fixed cells---')
            name = fl[fl.find('co/')+3:fl.find('.npy')]+'_changed_list.npy'
            newname = f'/media/HlabShare/models/inputs_new_donotdelete_updated/{name}'
            np.save(newname, tofix)
        
    elif keepgoing.lower() in ['n', 'no', 'nah', 'fuck you', 'fuck right off', 'i hate you', 'i hate myself', 'plz no']:
        print("\nYea I feel that. bye")
        sys.exit()
    else:
        print('\nk i dont know how to take that. so were just gonna do another file\n')