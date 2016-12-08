import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from definitions import ZERO_DIVISION_SAFE


class Incremental_SAX:
  """Class gathering all required function to run SAX during an online acquisition of time series
  We suppose that we already have a certain amount (frame_size) of points"""
    def __init__(self, frame_size , alphabet_size, nb_pieces):
        self.window_size = frame_size
        self.alphabet_size = alphabet_size
        self.nb_subdivision = nb_pieces
        self.global_mean = 0
        self.global_variance = 0 + ZERO_DIVISION_SAFE
        self.global_frequency = 0
        self.ts = []

    def update_window(self,new_point):
        


    