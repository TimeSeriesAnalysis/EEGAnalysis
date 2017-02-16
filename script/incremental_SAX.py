import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from definitions import ZERO_DIVISION_SAFE


# TO DO : CHECK IF WE HAVE A COLUMN ARRAY

class Incremental_SAX:
    """
        Class gathering all required function to run SAX during an online acquisition of time series.
        We assume that we already have a certain amount of points (500). The result of the algorithm is stored
        in the attribute called SAX.
        :param alphabet_size : Length of the alphabet which is numeric
        :type alphabet_size : Integer
        :param nb_subwindow : Number of sub-window (discretization) composing the sliding window
        :type nb_subwindow : Integer
        :param length_subwindow : Length of each sub-window(discretization) composing the sliding window
        :type length_subwindow : Integer
        :param time_serie : Data to transform
        :type time_serie : Numpy array of time series (floats numbers)
    """
    def __init__(self, alphabet_size, nb_subwindow , length_subwindow, time_serie):
        self.alphabet_size = alphabet_size
        self.alphabet = range(alphabet_size)
        self.nb_subwindow = nb_subwindow
        self.len_subwindow = length_subwindow
        self.window_size = self.nb_subwindow * self.len_subwindow                     #The sliding window size is defined as the product of the number of sub-window and its size
        self.oldest = 0                                                                                                #Define the first point which will be removed of the sliding frame (the oldest one) 
        self.global_mean = self.window.mean(axis = 0)                                          #Mean of each series contained in the sliding window
        self.global_variance = self.window.var(axis = 0) + ZERO_DIVISION_SAFE  #Variance of each series contained in the sliding window
        self.global_frequency = self.window_size                                                      #Frequency used to normalize the data each time we have a new point
        self.subwin_means = np.asarray(map(lambda xs: xs.mean(axis = 0), np.array_split(self.window, self.nb_subdivision)))         #Numpy array containing the mean of each subwindow 
        self.percentile = np.percentile((time_serie - time_serie.mean())/time_serie.std(),np.linspace(1./alphabet_size, 1-1./alphabet_size, alphabet_size-1)*100)       #percentiles of the normalized initial data
        self.znormalization()
        self.SAX = np.asarray([(self.alphabet[0] if ts_value < self.percentile[0] else (self.alphabet[-1] if ts_value > self.percentile[-1] else self.alphabet[np.where(self.percentile <= ts_value)[0][-1]+1])) for ts_value in self.window])          #Compute SAX transformation on initial data 
        self.unormalization()


    def update_window(self,new_point):
        """
            Add the new collected point and remove the oldest one in window attribute. It also updates the sliding window's features (mean, frequency and variance)
            :param new_point : New point collected to add
            :type new_point : Float number
        """
        self.window[self.oldest] = new_point
        new_global_frequency = self.global_frequency + 1
        new_global_mean = (self.global_mean * self.global_frequency + new_point) * 1./ (new_global_frequency)
        self.global_variance = (self.global_frequency * (self.global_variance + (self.global_mean - new_global_mean)**2) + (new_point - new_global_mean)**2) * 1./(new_global_frequency)
        self.global_mean = new_global_mean
        self.global_frequency = new_global_frequency


    # def update_global_mean_variance(self, new_point):                 Just a temporary solution
    #     new_global_frequency = self.global_frequency + 1
    #     new_global_mean = (self.global_mean * self.global_frequency + new_point) * 1./ (new_global_frequency)
    #     self.global_variance = (self.global_frequency * (self.global_variance + (self.global_mean - new_global_mean)**2) + (new_point - new_global_mean)**2) * 1./(new_global_frequency)
    #     self.global_mean = new_global_mean
    #     self.global_frequency = new_global_frequency


    def znormalization(self):
        self.window = (self.window - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)

    def PAA_SAX_transform(self):
        #local_subwin = np.copy(self.subwin_means)
        for index,mean in enumerate(self.subwin_means) :
            oldest_first = self.oldest + self.sublen * index
            oldest_next = self.oldest + self.sublen * (index + 1)
            if oldest_first >= self.window_size :
                oldest_first -= self.window_size
            if oldest_next >= self.window_size :
                oldest_next -= self.window_size
            local_subwin_mean = (mean*self.sublen - self.window[oldest_first] + self.window[oldest_next]) * 1.0 / self.sublen
            if local_subwin_mean > self.subwin_means[index] :
                if self.SAX[index] > self.percentile[-1]:
                    self.SAX[index] = self.percentile[-1]
                elif local_subwin_mean > self.percentile[self.SAX[index]] :
                    self.SAX[index] = self.alphabet[np.where(self.percentile[self.SAX[index]:] <= local_subwin_mean)[0][-1]]
            elif local_subwin_mean < self.subwin_means[index] :
                if self.SAX[index] - 1 < 0:
                    self.SAX[index] = 0
                elif local_subwin_mean < self.percentile[self.SAX[index] - 1]  :
                    self.SAX[index] = self.alphabet[np.where(self.percentile[:self.SAX[index]] <= local_subwin_mean)[0][-1]]
            self.subwin_means[index] = local_subwin_mean
        self.oldest += 1
        if self.oldest >= self.window_size : 
            self.oldest = 0

    def unormalization(self) :
        self.window =  self.window * np.sqrt(self.global_variance) + self.global_mean

    def run(self, new_point):
        self.update_window(new_point)
        self.update_global_mean_variance(new_point)
        self.znormalization()
        self.PAA_SAX_transform()
        self.unormalization()

if __name__ == "__main__":
    isax = Incremental_SAX(10,10,10,np.random.randn(500))
    print isax.percentile
    print isax.SAX
    print isax.window
    isax.run(10)
