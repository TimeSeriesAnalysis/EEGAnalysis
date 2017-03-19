# Import already existing modules
import numpy as np

# Import code from our own files
from definitions import ZERO_DIVISION_SAFE




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
    def __init__(self, alphabet, nb_subwindow, length_subwindow, time_series):
        self.alphabet = alphabet
        self.alphabet_size = length(alphabet)
        self.nb_subwindow = nb_subwindow
        self.len_subwindow = length_subwindow
        self.window_size = self.nb_subwindow * self.len_subwindow                     #The sliding window size is defined as the product of the number of sub-window and its size
        self.window = time_series[:,:self.window_size]                                                #Creation of the sliding window containing raw data. It must be a numpy array
        self.sorted_distribution = np.sort(self.window, axis = 1, kind = 'mergesort')
        self.index_oldest = 0                                                                                      #Define index of the first point to be removed when we get a new point 
        self.global_mean = self.window.mean(axis = 1)                                          #Mean of each series contained in the sliding window
        self.global_variance = self.window.var(axis = 1) + ZERO_DIVISION_SAFE  #Variance of each series contained in the sliding window
        self.percentils_index = [i * (self.window_size - 1) / self.alphabet_size for i in xrange(self.alphabet_size - 1)]
        self.znormalization()
        self.percentils = [self.window[i] for i in self.percentils_index]
        self.subwin_means = np.asarray(map(lambda xs: xs.mean(axis = 1), np.hsplit(self.window, self.nb_subwindow))).T         #Numpy array containing the mean of each subwindow 
        self.SAX = np.asarray([(self.alphabet[0] if ts_value < self.percentils[0] else (self.alphabet[-1] if ts_value > self.percentils[-1] else self.alphabet[np.where(self.percentils <= ts_value)[0][-1] + 1])) for ts_value in self.subwin_means])          #Compute a SAX transformation on initial data 
        self.unormalization()


    def update_window(self,new_point):
        """
            Add the new collected point and remove the oldest one in window attribute. It also updates the sliding window's features (mean, frequency and variance)
            :param new_point : New point collected to add
            :type new_point : Float number
        """
        removed_point = self.window[self.index_oldest]
        temp_mean = self.global_mean
        self.global_mean = temp_mean + (new_point - removed_point) * 1. / self.window_size
        self.global_variance = self.global_variance + (new_point**2 -removed_point**2 + 2*temp_mean*(new_point - removed_point) + (new_point - removed_point)**2 * 1. / self.window_size) *1. /self.window_size
        self.window[self.index_oldest] = new_point
        


    def znormalization(self):
        """
            Normalize the whole window by substracting its own mean and dividing by its own variance.
        """
        self.window = (self.window - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)
        self.sorted_distribution  = (self.sorted_distribution - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)


    def PAA_SAX_transform(self):            #To have good quantiles we'll need to find good index for each quntile ---> ordered list splitted in alpha parts then tkae the goods one to have quantiles ! Then to update just need to delete one and search where to put the oher
        """
            Perform the discretization and store the transformation into SAX attribute.
        """
        for index, mean in enumerate(self.subwin_means) :
            oldest_first = self.oldest + self.len_subwindow * index
            oldest_next = self.oldest + self.len_subwindow * (index + 1)
            if oldest_first >= self.window_size :
                oldest_first -= self.window_size
            if oldest_next >= self.window_size :
                oldest_next -= self.window_size
            local_subwin_mean = (mean*self.len_subwindow - self.window[oldest_first] + self.window[oldest_next]) * 1.0 / self.len_subwindow
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
        """
            Perform the exact opposite of znormalization function and un-normalize sliding frame's content.
        """
        self.window =  self.window * np.sqrt(self.global_variance) + self.global_mean
        self.sorted_distribution =  self.sorted_distribution * np.sqrt(self.global_variance) + self.global_mean


    def run(self, new_point):
        """
            Function gathering all required stage to perform a proper incremental SAX transformation.
            You need to call this function each time you have a new point.
            :param new_point : New collected point to add
            :type new_point : Float number
        """
        self.update_window(new_point)
        #self.update_global_mean_variance(new_point)        #useless since the merge of two functions
        self.znormalization()
        self.PAA_SAX_transform()
        self.unormalization()