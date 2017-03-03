# Objective

Our goal is to create or use algorithm enabling us to classify or to clusterize EEG signals for medical purposes (epilepsy detection ...).
The long term goal is to make our approaches more general and use them to analyse time series in different fields.

# [Symbolic Aggregate approXimation Algorithm (SAX)](https://timeseriesanalysis.github.io/static_SAX.html)
The aim of this algorithm is to convert a time serie into a sequence of letters. An incremental version of SAX has also been produced
If this approach is not sufficient we plan to use short time Fourier transform to map locally our time series
into a frequency representation and use more traditional methods

# [Classifying pieces of time series using Hidden Markov Model (HMM)](https://timeseriesanalysis.github.io/HMM.html)
Apply an HMM to the sequence of letters to produce a probabilistic model of the sequence of letters obtained from the EEG signal

# Classifying pieces of time series using Recurrent Neural Networks
Apply a RNN to the sequence of letters to produce a model of the sequence of letters obtained from the EEG signal

# Useful links
## DATASETS
+ [UCI eeg eyes](https://archive.ics.uci.edu/ml/datasets/EEG+Eye+State)

## Related topics
+ [neurofeeback](https://en.wikipedia.org/wiki/Neurofeedback)
+ [Event-related potential](https://en.wikipedia.org/wiki/Event-related_potential)

## Hidden Markov Models
+ [ghmm](http://ghmm.org)
+ [blog post](http://aimotion.blogspot.fr/2011/05/hidden-markov-models.html)

## Dynamic monitoring quantiles
+ [Dynamic monitoring quantiles stackoverflow](http://stats.stackexchange.com/questions/7959/algorithm-to-dynamically-monitor-quantiles)

## EEG
+ [EEG](http://emedicine.medscape.com/article/1139332-overview#a2)

# TODO
+ Implement Dynamic monitoring quantiles to the incremental version of SAX
+ Test if SAX transformation is sufficient or FFT transformation is required
+ Evaluate and compare HMM results and RNN 
+ Test on real time EEG signals

