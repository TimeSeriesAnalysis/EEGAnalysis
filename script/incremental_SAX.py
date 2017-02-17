# Import already existing modules
import numpy as np

# Import code from our own files
from definitions import ZERO_DIVISION_SAFE




class Incremental_SAX:
    """
        Class gathering all required function to run SAX during an online acquisition of time series.
        We assume that we already have a certain amount of points chosen by the user through two following parameters. The result of the algorithm is stored
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
    def __init__(self, alphabet_size, nb_subwindow, length_subwindow, time_serie):
        self.alphabet_size = alphabet_size
        self.alphabet = range(alphabet_size)
        self.nb_subwindow = nb_subwindow
        self.len_subwindow = length_subwindow
        self.window_size = self.nb_subwindow * self.len_subwindow                     #The sliding window size is defined as the product of the number of sub-window and its size
        self.window = np.asarray(time_serie[:self.window_size])                            #Creation of the sliding window containing raw data. It must be a numpy array
        self.oldest = 0                                                                                                #Define index of the first point removed when updating the sliding frame (the oldest one) 
        self.global_mean = self.window.mean(axis = 0)                                          #Mean of each series contained in the sliding window
        self.global_variance = self.window.var(axis = 0) + ZERO_DIVISION_SAFE  #Variance of each series contained in the sliding window
        self.global_frequency = self.window_size                                                      #Frequency used to normalize the data each time we have a new point
        self.subwin_means = np.asarray(map(lambda xs: xs.mean(axis = 0), np.array_split(self.window, self.nb_subwindow)))         #Numpy array containing the mean of each subwindow 
        self.percentile = np.percentile((self.window - self.global_mean) / self.global_variance, np.linspace(1. / self.alphabet_size, 1 - 1. / alphabet_size, alphabet_size - 1) * 100)       #percentiles of the normalized initial data
        # self.znormalization()
        # self.SAX = np.asarray([(self.alphabet[0] if ts_value < self.percentile[0] else (self.alphabet[-1] if ts_value > self.percentile[-1] else self.alphabet[np.where(self.percentile <= ts_value)[0][-1] + 1])) for ts_value in self.window])          #Compute SAX transformation on initial data 
        # self.unormalization()


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
        """
            Normalize the whole window by substracting its own mean and dividing by its own variance.
        """
        self.window = (self.window - self.global_mean) / (np.sqrt(self.global_variance)+ZERO_DIVISION_SAFE)


    def PAA_SAX_transform(self):
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


    def run(self, new_point):
        """
            Function gathering all required stage to perform a proper incremental SAX transformation.
            You need to call this function each time you have a new point.
            :param new_point : New collected point to add
            :type new_point : Float number
        """
        #self.update_window(new_point)
        #self.update_global_mean_variance(new_point)        #useless since the merge of two functions
        self.znormalization()
        self.PAA_SAX_transform()
        self.unormalization()


    def see(self):
        """
            Enabling to visualize step by step
        """
        #self.update_window(new_point)
        #self.update_global_mean_variance(new_point)        #useless since the merge of two functions
        self.znormalization()
        self.PAA_SAX_transform()
        self.unormalization()

if __name__ == "__main__":
    isax = Incremental_SAX(10,10,10,np.random.randn(500))
    print "INITIAL DATA : ", isax.window
    isax.see()
    print "SAX TRANSFORMATION : ", isax.SAX
    # print isax.percentile
    # print isax.SAX
    # print isax.window
    #isax.run(np.random.randn)
