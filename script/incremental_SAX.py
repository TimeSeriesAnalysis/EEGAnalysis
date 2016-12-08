import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from definitions import ZERO_DIVISION_SAFE


class Incremental_SAX:
  """Class gathering all required function to run SAX during an online acquisition of time series
  We assume that we already have a certain amount (frame_size) of points"""
    def __init__(self, alphabet_size, nb_pieces ,nb_inside_segment,time_serie):
        self.alphabet_size = alphabet_size
        self.nb_subdivision = nb_pieces
        self.window_size = nb_pieces * nb_inside_segment
        self.ts = np.asarray(time_serie)
        self.global_mean = ts.mean(axis = 0)
        self.global_variance = ts.var(axis = 0) + ZERO_DIVISION_SAFE
        self.global_frequency = self.window_size
        
    def update_global_mean_variance(self, x):
        new_global_frequency = self.global_frequency + 1
        new_global_mean = (self.global_mean * self.global_frequency + x) * 1./ (new_global_frequency)
        self.global_variance = (self.global_frequency * (self.global_variance + (self.global_mean - new_global_mean)**2) + (x - new_global_mean)**2) * 1./(new_global_frequency)
        self.global_mean = new_global_mean
        self.global_frequency = new_global_frequency


        


    