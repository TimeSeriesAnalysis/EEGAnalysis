import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from definitions import ZERO_DIVISION_SAFE

def znormalization(ts):
    mus = ts.mean(axis = 0)
    stds = ts.std(axis = 0)
    return (ts - mus) / (stds+ZERO_DIVISION_SAFE)

def paa_transform(ts, n_pieces):
    """
    ts: the columns of which are time series represented by e.g. np.array
    n_pieces: M equally sized piecies into which the original ts is splitted
    """
    splitted = np.array_split(ts, n_pieces) ## along columns as we want
    return np.asarray(map(lambda xs: xs.mean(axis = 0), splitted))

def sax_transform(ts, n_pieces, alphabet_sz):
    """
    ts: columns of which are time serieses represented by np.array
    n_pieces: number of segments in paa transformation
    alphabet: the letters to be translated to, e.g. "abcd", "ab"
    return np.array of ts's sax transformation
    Steps:
    1. znormalize
    2. ppa
    3. find norm distribution breakpoints by scipy.stats
    4. convert ppa transformation into strings
    """
    alphabet = range(alphabet_sz)
    def translate(ts_values):
        return np.asarray([(alphabet[0] if ts_value < thrholds[0]
                else (alphabet[-1] if ts_value > thrholds[-1]
                      else alphabet[np.where(thrholds <= ts_value)[0][-1]+1]))
                           for ts_value in ts_values])
        
    ts_norm = znormalization(ts)
    thrholds = np.percentile(ts_norm,np.linspace(1./alphabet_sz, 
                                          1-1./alphabet_sz, 
                                          alphabet_sz-1)*100)
    paa_ts = paa_transform(ts_norm, n_pieces)
    return np.apply_along_axis(translate, 0, paa_ts)