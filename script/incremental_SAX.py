import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from definitions import ZERO_DIVISION_SAFE


# TO DO : CHECK IF WE HAVE A COLUMN ARRAY

class Incremental_SAX:
  """Class gathering all required function to run SAX during an online acquisition of time series
  We assume that we already have a certain amount (frame_size) of points"""
    def __init__(self, alphabet_size, nb_pieces ,nb_inside_segment,time_serie):
        self.alphabet_size = alphabet_size
        self.nb_subdivision = nb_pieces
        self.sublen = nb_inside_segment
        self.window_size = nb_pieces * self.sublen
        self.window = np.asarray(time_serie)
        self.oldest = self.window[0]
        self.global_mean = ts.mean(axis = 0)
        self.global_variance = ts.var(axis = 0) + ZERO_DIVISION_SAFE
        self.global_frequency = self.window_size
        self.subwin_means = np.asarray(map(lambda xs: xs.mean(axis = 0), np.array_split(self.window, self.nb_subdivision)))

    def update_window(self,new_point):
        self.window[self.oldest] = new_point
      
    def update_global_mean_variance(self, x):
        new_global_frequency = self.global_frequency + 1
        new_global_mean = (self.global_mean * self.global_frequency + x) * 1./ (new_global_frequency)
        self.global_variance = (self.global_frequency * (self.global_variance + (self.global_mean - new_global_mean)**2) + (x - new_global_mean)**2) * 1./(new_global_frequency)
        self.global_mean = new_global_mean
        self.global_frequency = new_global_frequency

    def znormalization(self):
        return (self.window - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)

    def paa_transform(self):
        for index,mean in enumerate(self.subwin_means) :
            oldest_first = self.oldest + self.sublen * index
            oldest_next = self.oldest + self.sublen * (index + 1)
            if oldest_first > self.window_size :
                oldest_first -= self.window_size
            if oldest_next > self.window_size :
                oldest_next -= self.window_size
            self.subwin_means[index] = (mean*self.sublen - self.window[oldest_first] + self.window[oldest_next + self.sublen]) * 1.0 / self.sublen
        self.oldest += 1
        if self.oldest >= self.window_size : 
            self.oldest = 0

    def unormalization(self) :
        return self.ts * np.sqrt(self.global_variance) + self.global_mean


        


    