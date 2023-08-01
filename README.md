# criticality
This repository contains the code, written in Matlab, that computes avalanches, and plots the histogram of avalanches' sizes and durations, both with linear and logarithmic binning 
It refers to  paper  "Criticality of neuronal avalanches in human sleep and their relationship with sleep macro- and micro-architecture".

The code ok_avalanches.m computes the avalanches. 
An avalanche is defined as a continuous time interval in which there is at least one amplitude excursion beyond threshold in at least one EEG channel. Avalanches are preceded and followed by time intervals with no excursions beyond threshold on any EEG channel.
In the code artifact-free EEG signals are z-score normalized to have zero mean and unit standard deviation. To define the threshold θ, we analyzed the distribution of EEG amplitudes, to check when they deviate from the gaussian. For each EEG channel, large positive or negative excursions beyond the threshold θ = ± n SD were identified. The size of an avalanche, s, was defined as the sum over all channels of the absolute values of the signals exceeding the threshold.
The code  gives 
1) a figure of the probability density of the z-score normalized EEG signal amplitude, with fitted gaussian.
2) the histogram of avalanche's size and duration distributions, both with linear and logarithmic binning.
The code ok_avalanches.m load data from file SS_N1_TEST.mat that have a test sample of signals (only the few hours of the eeg signals of first subject, where waking and motion artifact segments were removed).
   
    Prerequisites: 
        Matlab
