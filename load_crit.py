import numpy as np 
from sahara_work import Crit
import matplotlib.pyplot as plt 

path = ''
crit = np.load(path, allow_pickle=True)[0]

# stats are all stored in the object like this
dcc = crit.dcc
p_value_burst = crit.p_value_burst

#figures are also stored
scaling_fig = crit.scaling_plot

#to view the figure 
crit.scaling_plot.show()
