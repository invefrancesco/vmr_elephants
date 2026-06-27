This repository provides an R implementation for the estimation and simulation of Circular Mixture Autoregressive models based on the von Mises distribution. The model is specifically designed to handle circular time-series data across multiple mixture components and independent data bursts.

The files are organized as follows: 

* `fun_estimation.R` contains the log-likelihood function, the Q function and th EM algorithm
* `fun_simulation.R` contains the functions for generating a sigle burst of angles.
* `main.R` run a full workflow, including parameter setup, data simulation, model fitting, and results visualization
* `charts.R`contains the visualization suite for the project

#### Q
[ ] "Brucio" 500 osservazioni restando nello stesso stato e poi genero le ultime dallo stato estratto con multinomial logit (?);