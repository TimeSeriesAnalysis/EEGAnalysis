# Import already existing modules
import numpy as np
from bisect import insort_left as sinsert
from bisect import bisect_left as sindex
# Import code from our own files
from definitions import *




class Dynamic_SAX:
    """
        Class gathering all required function to run SAX during an online acquisition of time series.
        Everything is implemented to handle multi-dimensional signals.
        We assume that we already have a certain amount of points chosen by the user through two following parameters. The result of the algorithm is stored in the attribute called SAX.
        :param alphabet : Desired alphabet used to transform your signal
        :type alphabet : List of whatever you want
        :param nb_subwindow : Number of sub-window (interval of discretization) composing the sliding window
        :type nb_subwindow : Integer
        :param length_subwindow : Length of each sub-window(interval of discretization) composing the sliding window
        :type length_subwindow : Integer
        :param time_series : Data to transform
        :type time_series : Numpy array of time series (floats numbers) each row must be a time serie.
    """
    def __init__(self, alphabet, nb_subwindow, length_subwindow, time_series, option = VERTICAL):
        self.alphabet = alphabet
        self.alphabet_size = len(alphabet)
        self.nb_subwindow = nb_subwindow
        self.len_subwindow = length_subwindow
        self.window_size = self.nb_subwindow * self.len_subwindow                     #The sliding window size is defined as the product of the number of sub-window and its size
        self.window = time_series[:,:self.window_size].T                                                #Creation of the sliding window containing raw data. It must be a numpy array
        self.sorted_distribution = np.sort(self.window, axis = 0, kind = 'mergesort')
        self.index_oldest = 0                                                                                      #Define index of the first point to be removed when we get a new point 
        self.global_mean = self.window.mean(axis = 0)                                          #Mean of each series contained in the sliding window
        self.global_variance = self.window.var(axis = 0) + ZERO_DIVISION_SAFE  #Variance of each series contained in the sliding window
        self.percentils_index = np.asarray([i * (self.window_size - 1) / self.alphabet_size for i in xrange(self.alphabet_size - 1)])
        self.znormalization()
        self.percentils = self.sorted_distribution[self.percentils_index,:]#[self.sorted_distribution[i,:] for i in self.percentils_index]
        self.subwin_means = np.asarray(map(lambda xs: xs.mean(axis = 0), np.vsplit(self.window, self.nb_subwindow)))         #Numpy array containing the mean of each subwindow 
        self.SAX = np.zeros(self.subwin_means.shape)
        self.sorted_distribution = self.sorted_distribution.T.tolist()
        self.percentils = self.percentils.T.tolist()
        self.PAA_to_SAX()
        #self.unormalization()


    def update_global_parameters(self,new_point):
        """
            Add the new collected point and remove the oldest one in window attribute. It also updates the sliding window's features (mean, frequency and variance)
            :param new_point : New point collected to add
            :type new_point : Float number
        """
        new_point = (new_point - self.global_mean) * 1.0 / self.global_variance 
        removed_point = self.window[self.index_oldest]
        temp_mean = self.global_mean
        self.global_mean = temp_mean + (new_point - removed_point) * 1. / self.window_size
        temp_var =  self.global_variance
        self.global_variance = temp_var + (new_point**2 -removed_point**2 - 2*temp_mean*(new_point - removed_point) - (new_point - removed_point)**2 * 1. / self.window_size) *1. /self.window_size        
        for index,distribution in enumerate(self.sorted_distribution):
            self.sorted_distribution[index].remove(removed_point[index])
            sinsert(self.sorted_distribution[index], new_point[index])
        self.update_percentils()

    def update_percentils(self):
        for i,percentil in enumerate(self.percentils):
            for j,value in enumerate(percentil):
                self.percentils[i][j] = self.sorted_distribution[i][self.percentils_index[j]] 

    def znormalization(self):
        """
            Normalize the whole window by substracting its own mean and dividing by its own variance.
        """
        self.window = (self.window - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)
        self.sorted_distribution  = (self.sorted_distribution - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)


    def PAA_transform(self, new_point):            #To have good quantiles we'll need to find good index for each quntile ---> ordered list splitted in alpha parts then tkae the goods one to have quantiles ! Then to update just need to delete one and search where to put the oher
        """
            Perform the discretization and store the transformation into SAX attribute.
        """
        for index, mean in enumerate(self.subwin_means):
            current_oldest = self.index_oldest + self.len_subwindow * index
            if current_oldest >= self.window_size :
                current_oldest -= self.window_size
            self.subwin_means[index] -= self.window[current_oldest] * 1./self.len_subwindow
        self.window[self.index_oldest] = (new_point - self.global_mean) * 1.0 / self.global_variance
        for index, mean in enumerate(self.subwin_means):
            next_oldest = self.index_oldest + self.len_subwindow * (index+1)
            if next_oldest >= self.window_size :
                next_oldest -= self.window_size
            self.subwin_means[index] -= self.window[next_oldest] * 1./self.len_subwindow
        self.index_oldest += 1
        if self.index_oldest == self.window_size : 
            self.index_oldest = 0

    def PAA_to_SAX(self):
        for i in xrange(self.SAX.shape[0]):
            for j in xrange(self.SAX.shape[1]):
                self.SAX[i][j] = self.alphabet[sindex(self.percentils[j],self.subwin_means[i][j])]

    def SAX_transform(self, new_point):
        self.update_global_parameters(new_point)
        self.PAA_transform(new_point)
        self.PAA_to_SAX()
