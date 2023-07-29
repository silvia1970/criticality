# criticality
This repository contains the code, written in Matlab, that computes avalanches, and plots the histogram of avalanches' sizes and durations, both with linear and logarithmic binning 
It refers to  paper  "Criticality of neuronal avalanches in human sleep and their relationship with sleep macro- and micro-architecture".

The code ok_N1_avalanches.m  loads data from 'SS_N1.mat' and computes avalanches.
The file SS_N1.mat contains the EEG signals of first subject (where waking and motion artifact segments were removed).
The code ok_N1_avalanches.m gives 
1) a figure of the probability density of the z-score normalized EEG signal amplitude, and 
2) the histogram of avalanche's size and duration distributions, both with linear and logarithmic binning.
   
    Prerequisites: 
        Matlab
