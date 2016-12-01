# starting

First step of test and exploration of methods to analyse different kind of curves.
Our goal is to create or use algorithm enabling us to classify or to clusterize EEG for medical purposes. But if we succeed, the goal would be to get that more general and use that kind of approch in different fields.

# IDEAS 
## First operations applied to the time series
+ Use of [SAX](https://github.com/dolaameng/pysax) in order to convert our time series into a sequence of letters.
+ If this approach is not sufficient we plan to use ?Short time Fourier transform to map locally our time series into a frequency representation and use more traditional methods

## Classifying pieces of time series  using HMM
+ Apply an HMM to the sequence of letters. 

## STDP neural network
+ See with Ilya if such a method could be applied


# TODO

+ Search more information about home made / cheap EEG
+ Possible applications
+ Evaluate and compare different approaches
+ Exploring and coding 1d-SAX to compare the two algorithm
+ Code our HMM to classify data depending on pattern
+ Start building "online" version

#UNFINISHED TASKS
+ Code our own SAX generalist method 

# TO EXPLORE
+ [Time series Analysis](https://sflscientific.com/data-science/) 
+ [Kaggle challenge Winner on EEG prediction](https://www.kaggle.com/c/grasp-and-lift-eeg-detection/)
+ [ARMA](https://bicorner.com/2015/11/16/time-series-analysis-using-ipython/)

# TOOLS ?
+ [Riemann](https://github.com/alexandrebarachant/pyRiemann)
+ [Dynamic time wrapping](https://en.wikipedia.org/wiki/Dynamic_time_warping)
+ [Python library for EEG Analysis](http://ptsa.readthedocs.io/en/latest/index.html)
