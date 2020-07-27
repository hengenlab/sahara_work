import numpy as np
from copy import deepcopy as cdc

def get_avalanches(data, perc = 0.25):
    '''
        Function that goes through an array of binned spikes and determines the avalanch boundaries and properties

        parameters:
            data - array of spike times. one row for each neuron
            perc - threshold for defining an avalanche, if network is dominated by silent periods you can use 0

        returns:
            Result - a dictionary with 2 inputs. 'S' is the size of each avalanche (number of spikes above threshold).
                    'T' is the duration (number of time bins avalanche spanned)
    '''
    n, m = np.shape(data) # num cells, num bins

    if n == 1:
        network = cdc(data)
    else:
        network = np.nansum(data, axis = 0) # collapse into single array. sum the amount of activity in each bin.

    if perc > 0:
        sortN = np.sort(network)
        threshold2 = sortN[round(m * perc)] # determine the treshold. if perc is .25, then its 25% of network activity essentially
    else:
        threshold2 = 0

    zdata = cdc(network)
    z2data = cdc(network)

    zdata[zdata <= threshold2] = 0 # intervals
    zdata[zdata > threshold2] = 1 # avalanches

    Z = np.where(zdata == 0)[0] # location of intervals

    z1data = np.delete(zdata, Z[np.where(np.diff(Z) == 1)[0]])  # cuts out irrelevant 0s, series of 1s and a single 0 separation
    avalanches = np.delete(z2data, Z[np.where(np.diff(Z) == 1)[0]])
    avalanches[z1data == 0] = 0  # use a single 0 to separate network activities in each avalanche

    J = np.where(z1data == 0)[0]  # location of the intervals

    burst = []
    for i in np.arange(0, np.size(J) - 1):
        # print(J[i+1]-J[i]-2)
        fired = np.sum(avalanches[J[i] + 1:J[i + 1]]) - threshold2 * (J[i + 1] - J[i] - 1)
        burst.append(fired)

    T = np.diff(J) - 1  # AVduration
    T = T[T > 0]  # Duration should be positive

    Result = {
        'S': np.asarray(burst),
        'T': T
    }

    return Result
