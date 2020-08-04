import numpy as np
import neuraltoolkit as ntk 

trode_pairs = []
for i in np.arange(0,64):
    chdat = d[i]
    corrcoefs = [np.corrcoef(chdat, d[x])[0][1] for x in np.arange(0,64)]
    s = np.sort(corrcoefs)
    idxs = np.argsort(corrcoefs)
    trode = np.sort(idxs[-4:]) # sort goes smallest to largest, so we want the 4 largest numbers and since its linear, just the indexes works
    trode_pairs.append(list(trode))