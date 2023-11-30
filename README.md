# Diagnosis-of-Parkinson-s-Disease-using-Voice-Analysis

This project was done as part of my B.Sc Final Project with my teammate Shir Coffman.

The objective of this project is to use machine learning algorithms to identify Parkinson diseased based on changes in speech acoustic characteristics. 
The data set is based on 48 speakers, healthy and diseased subjects, each has approximately 10 different records: records of /AH/, /Pa-Ta-Ka/ and continuous speech.

Several previous studies have shown the ability to identify Parkinson's disease by analyzing the acoustic changes of speech and by using machine learning methods. 
Overall, this project extracts classical features that were noted in literature review as: HNR, NHR, SNR, Jitter and Shimmer, 
but also novel and original parameters that do not appear in previous studies. 
These parameters aim to evaluate the periodic nature of the signal and to quantify the level of the irregular frames, and with it, to identify "Parkinson voice". 
Such as: Special Frames, Autocorrelation Amplitudes feature, Entropy, 'S' features, Pitch and RMS decay, Spectral Slope, Speech Rate and Lengths Statistics

As part of the feature extraction process: 
1. Pre-processing of the records is done using MATLAB software, including division of the audio signal into frames and application of voice activity detection algorithm.
2. The acoustic characteristics are independently implemented, some based on previous articles, as jitter, shimmer, HNR and pitch. This phase includes statistical analysis of the features, if the feature returns value per frame.
3. Their separation capability between healthy and Parkinson's disease subjects is examined using machine learning classifiers of KNN and SVM while only those that provide high identification results are selected.

The final reesults and entire process is depicted in the Project Book.
