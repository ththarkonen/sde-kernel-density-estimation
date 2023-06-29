
clc
clear
close all

includeFolders = genpath('include');
addpath( includeFolders );

data = load("data/heston/heston.mat");

S = data.S;
V = data.V;

S = cumsum([log(97.0677); S]);
S = exp(S);

nData = length( S );

settings.dt = 0.003968253968254;
settings.simulationSteps = 1;
settings.r = 0.04;
settings.q = 0.015;
settings.nParticles = 15000;

dataHeston = struct();
dataHeston.S = S;
dataHeston.V = V;

ssFun = @(theta, data) -2 * hestonParticleLikelihood( data, theta, settings);

model = struct();
model.ssfun  = ssFun;
model.sigma2 = 1.0;

model.N = nData;

nIterations = 50000;

optionsAM.nsimu = nIterations;
optionsAM.updatesigma = 0;
optionsAM.waitbar = 1;
optionsAM.verbosity = 2;
optionsAM.method = 'dram';
optionsAM.burnintime = 300;
optionsAM.adaptint = 300;
optionsAM.drscale = 5;

parameters = {
    {'theta', 0.1}
    {'kappa', 3}
    {'xi', 0.25, 0}
    {'rho', -0.8, -1, 1}
    {'lambda1', 4.0}
};

tic
[~, chainAM] = mcmcrun( model, dataHeston, parameters, optionsAM);
toc
saveResults( "heston", chainAM, [S, V]);
