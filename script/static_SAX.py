# Import already existing modules
import numpy as np
from scipy.stats import norm

# Import our own code
from definitions import ZERO_DIVISION_SAFE


def znormalization(data):
    """
        Return normalized data (assumed to be a numpy array) by substracting the mean value and divinding by the standard deviation
        :param data : Time-series data. It can contain different signals (one per row)
        :type data : Numpy array of floats
        :returns: Normalized data
        :rtype: Numpy floats array of exactly data's dimension
    """
    mus = data.mean(axis = 1)
    stds = data.std(axis = 1)
    return (data - mus) / (stds+ZERO_DIVISION_SAFE)


def paa_transform(data, nb_interval):
    """
        Perform and return PAA transformation on data for an amount of interval.
        :param data : Time-series data. It can contain different signals (one per row)
        :type data : Numpy array of floats
        :param nb_interval : Number of discretization over the whole data.
        :type nb_interval : Integer
        :returns : PAA transform (mean value of each interval)
        :rtype : Numpy array of floats
    """
    splitted = np.array_split(data, nb_interval) ## along columns as we want
    return np.asarray(map(lambda xs: xs.mean(axis = 0), splitted))


def sax_transform(ts, n_pieces, alphabet_sz, use_gaussian_assuption = False):
    """
    ts: columns of which are time serieses represented by np.array
    n_pieces: number of segments in paa transformation
    alphabet_size : number of segment on vertical axis (size of alphabet)
    Steps:
    1. znormalize
    2. ppa
    3. find norm distribution breakpoints by scipy.stats
    4. convert ppa transformation into strings
    """
    alphabet = range(alphabet_sz) # we choose here a numeric alphabet
    def translate(ts_values):
        return np.asarray([(alphabet[0] if ts_value < thrholds[0]
                else (alphabet[-1] if ts_value > thrholds[-1]
                      else alphabet[np.where(thrholds <= ts_value)[0][-1]+1]))
                           for ts_value in ts_values])

    ts_norm = znormalization(ts)
    quantiles = np.linspace(1./alphabet_sz, 1-1./alphabet_sz, alphabet_sz-1)
    if use_gaussian_assuption:
        thrholds = norm.ppf(quantiles)
    else:
        thrholds = np.percentile(ts_norm,quantiles*100)
    paa_ts = paa_transform(ts_norm, n_pieces)
    return np.apply_along_axis(translate, 0, paa_ts)