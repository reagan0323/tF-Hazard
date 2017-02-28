# tF-Hazard
Two-way functional hazards model


The package provides several useful tools for two-way functional data analysis. 
The tFSVD function performs two-way functional singular value decomposition. It decomposes a data matrix, which can be viewed as the discretization of some two-way smooth function with noise, into low-rank smooth components.

The tFHazards function estimates a two-way smooth hazard function as a function of two timescales (e.g., birthtime and lifetime, or time to event and calendar time). It takes in censored time-to-event data with birth time information, and outputs a smooth and low-rank two-way hazard function estimation.

The example.m file gives a survival data example to illustrate the usage of different methods.

Both tFSVD and tFHazards are useful exploratory tools in functional data analysis. tFSVD is generally applicable to two-way functional data such as mortality rate analysis (i.e., age vs year); tFHazards is applicable to two-way censored data such as call center waiting time analysis (i.e., arrival time vs waiting duration).

