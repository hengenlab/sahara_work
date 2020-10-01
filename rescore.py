import numpy as np
import musclebeachtools_hlab.musclebeachtools as mbt

def rescore(old, new):
    old_cells = np.load(old, allow_pickle = True)
    new_cells = np.load(new, allow_pickle = True)

    for i, cell in enumerate(old_cells):
        qual = cell.quality
        new_cells[i].set_qual(qual)

    np.save(old, new_cells)