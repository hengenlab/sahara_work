'''
%RLDECODE Run-length decoding of run-length encode data.
%
%   X = RLDECODE(LEN, VAL) returns a vector XLEN with the length of each run
%   and a vector VAL with the corresponding values.  LEN and VAL must have the
%   same lengths.
%
%   Example: rldecode([ 2 3 1 2 4 ], [ 6 4 5 8 7 ]) will return
%
%      x = [ 6 6 4 4 4 5 8 8 7 7 7 7 ];
%
%   See also RLENCODE.

%   Author:      Peter John Acklam
%   Time-stamp:  2002-03-03 13:50:38 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

---------------------
translated to python by sahara ensley 07/02/2020
'''
import numpy as np
def rldecode(l, val):
    # keep only runs whose length is positive
    k = l > 0
    l = l[k]
    val = val[k]

    # now perform the actual run-length decoding 
    i = np.cumsum(l, dtype=int)
    j = np.zeros(int(i[-1]))
    j[ i[0:-2]+1 ] = 1
    j[0] = 1
    x = val[np.cumsum(j, dtype=int)]

    return x