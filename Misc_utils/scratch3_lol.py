import glob 
import numpy as np 
import sahara_work as saw   
import musclebeachtools as mbt
import sys


files = glob.glob('/media/HlabShare/Clustering_Data/*/*/*/*/co/*neurons_group0.npy') 
files = [f for f in files if 'files_block' not in f]
animals = [saw.get_info_from_path(f)[0] for f in files]
animals = np.unique(animals)

toscore = []
for a in animals:
    geno = saw.get_genotype(a)
    if geno in ['te4', 'wt', 'e4']:
        up = a[:3].upper()+'000'+a[3:]
        probe = saw.get_probe(a, 'CA1')
        if probe == -1:
            print(a, ' no ca1')
        else:
            subset = sorted([f for f in files if up in f and probe in f])
            num_files = len(subset)
            beg = int(num_files * 0.1)
            end = int(num_files * 0.9)
            print(a, up, 'total: ', len(subset), ' first: ', beg, ' second: ', end)
            toscore.append(subset[beg])
            toscore.append(subset[end])


model = '/media/HlabShare/models/caf/xgb_model_k10nfrpr_prob2_0_corr10_caf_sah_basedev7'

for i, f in enumerate(toscore):
    print(f'DONE {i} FILES OF {len(toscore)}')
    print(f'------ STARTING {f}')

    keepgoing = input(f"\nWanna keep going or nah?  ")

    while keepgoing.lower() not in ['y', 'yes', 'sure', 'might as well', 'fuck it', 'ugh', 'n', 'no', 'nah', 'fuck you', 'fuck right off', 'i hate you', 'i hate myself', 'plz no']:
        keepgoing = input('\nok i know youre upset but lets keep it together shall we? wanna do another file? ')

    if keepgoing.lower() in ['y', 'yes', 'sure', 'might as well', 'fuck it', 'ugh']:
        tofix=[]
        cells = np.load(f, allow_pickle=True)
        #mbt.autoqual(cells,model)
        for cell in cells:
            if cell.quality < 4:
                ogq = cell.quality
                cell.checkqual()
                newq = cell.quality
                # if ogq != newq:
                #     print('\nNEW QUAL DIFFERENT - saving cell')
                #     tofix.append(cell)

        # if len(tofix) > 0:
        #     print('---saving fixed cells---')
        #     name = f[f.find('co/')+3:f.find('.npy')]+'_changed_list.npy'
        #     newname = f'/media/HlabShare/models/inputs_new_donotdelete_updated/{name}'
        #     np.save(newname, tofix)
        
        save = input('do you want to save this block for handscoring?')
        if save == 'y':
            newname = f[:f.find('.npy')] + '_HANDSCORED.npy'
            np.save(newname, cells)
    else:
        print('yup. about right.')
        sys.exit()