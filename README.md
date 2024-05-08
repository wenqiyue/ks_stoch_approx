# ks_stoch
Stochastic approximation of the Kuramoto-Sakaguchi model of coupled oscillator system


This repository contains the MATLAB scripts for numerically investigating the Kuramoto-Sakaguchi model of
coupled oscillators and performing a stochastic approximation whereby the effects of rogue oscillators on
synchronized oscillators is approximated with a 2-dim OU process. The parameters of the 2-dim OU process 
are obtained as best fits against data from the full original system. 

ks_sim_1.m: simulate an instance of the system for a short period of time. Perform some analysis

ks_sim_Delta.m: simulate the system over an extended period of time. Perform some analysis. Gather some data.
Obtaining long time trajectories. The long time trajectories are used for fitting and comparison.

plotting_ks_sim_1_Delta.m: Plot results from the long time simulation. Also perform fitting

ks_sim
