# Import already existing modules
import numpy as np
from scipy.stats import norm

# Import our own code
from definitions import ZERO_DIVISION_SAFE



def clean_row_outliers(row, ampl):
    mean = row.mean()
    sd = row.std()
    row[row < mean - 1.*ampl * sd] = mean - 1.*ampl * sd
    row[row > mean + 1.*ampl * sd] = mean + 1.*ampl * sd
    return row


def clean_outliers(data, ampl):
    """
        Change all values of data not includes in the intervale [mean - amplitude*sd ; mean +amplitude*sd] by the limits.
        :param data : Contains all the times series (one by row)
        :type data : Numpy array of floats
        :param ampl : Regulate the amplitude of the interval
        :type ampl : Float
        :returns : Filtered data
        :rtype : Numpy array of floats
    """
    return np.apply_along_axis(clean_row_outliers, 1, data, ampl)


def znormalization(data):
    """
        Return normalized data (assumed to be a numpy array) by substracting the mean value and divinding by the standard deviation
        :param data : Time-series data. It can contain different signals (one per row)
        :type data : Numpy array of floats
        :returns: Normalized data
        :rtype: Numpy floats array of exactly data's dimension
    """
    mus = data.mean(axis = 1).reshape([data.shape[0], 1])
    stds = data.std(axis = 1).reshape([data.shape[0], 1])
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
    splitted = np.hsplit(data, nb_interval)
    return np.asarray(map(lambda xs: xs.mean(axis = 1), splitted)).T



def paa_to_alphabet(ts_values, alphabet, quantiles):
    """
        Translate PAA transformed data into an alphabet.
        :param ts_values : PAA transformed data.
        :type ts_values : Numpy array of integers
        :param alphabet : Alphabet you want (can be numeric, alphanumeric or only letters)
        :type alphabet : List of what you wan
        :param quantiles : Quantiles of original data you want to use to translate.
        :type quantiles : List or numpy array of Floats
        :returns : Equivalent alphabet of input data
        :rtype : Numpy array of alphabet.
    """
    return np.asarray([(alphabet[0] if ts_value < quantiles[0]
            else (alphabet[-1] if ts_value > quantiles[-1]
                  else alphabet[np.where(quantiles <= ts_value)[0][-1]+1]))
                       for ts_value in ts_values])



def sax_transform(ts, n_pieces, alphabet_sz, use_gaussian_assumption = False):
    """
        Perform SAX transformation. First it normalizes data, then it applies PAA transformation and translate it into an alphabet (numeric here).
        :param ts : Times series to transform
        :type ts : Numpy array. Each row is a time serie.
        :param n_pieces : Number of segment in paa transformation
        :type n_pieces : Integer
        :param alphabet_sz : Length of the alhpabet you want (number of segment on the vertical axis)
        :type alphabet_sz : Integer
        :param use_gaussian_assumption : Tells if you want to use gaussian percentile or take real distribution percentile.
        :type use_gaussian_assumption : Boolean
        :returns : SAX transformed data
        :rtype : Numpy array of integers.
    """
    alphabet = range(alphabet_sz) # we choose here a numeric alphabet
    ts_norm = znormalization(ts)
    quantiles = np.linspace(1./alphabet_sz, 1-1./alphabet_sz, alphabet_sz-1)
    if use_gaussian_assumption:
        thrholds = norm.ppf(quantiles)
    else:
        thrholds = np.percentile(ts_norm,quantiles*100)
    paa_ts = paa_transform(ts_norm, n_pieces)
    return np.apply_along_axis(paa_to_alphabet, 1, paa_ts, (alphabet, thrholds))