# Kernel density estimation with Markov chain Monte Carlo (MCMC) for stochastic differential equations (SDE)
Implementation of MCMC sampling of SDE model parameters. This is a repository for the software used to create the kernel density estimate results for the paper 
"Bayesian synthetic likelihood for stochastic models with applications in mathematical finance'' ([doi.org/10.3389/fams.2023.1187878](https://doi.org/10.3389/fams.2023.1187878)).

The software approximates transition probality densities for Ornstein-Uhlenbeck, Merton, and Heston stochastic differential equation systems. The approximations are constructed by simulating SDE
trajectories for each of the systems with given parameters. The simulated trajectories are samples from the transition probality distributions. A kernel density estimate is fitted to the simulated trajectories
which is then used to compute an approximate likelihood for the model parameters. Finally, the approximate likelihood is used in MCMC sampling which yields samples from the posterior distributions of the SDE 
model parameters.

# Installation
Clone or download the repository and add the include folder and subfolders to your MATLAB path. This should happen automatically on running any of the scripts in the project root folder.

# References
If you find the software useful, please cite [doi.org/10.3389/fams.2023.1187878](https://doi.org/10.3389/fams.2023.1187878).
