import glob
import sahara_work as s 
from sahara_work import Crit

paths = glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/*/*/*/Crit*')

s.construct_fr_df(paths)