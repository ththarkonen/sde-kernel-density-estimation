
clc
clear
close all

includeFolders = genpath('include');
addpath( includeFolders );

data = load("data/OU/ou.mat");

x = data.X( 1, :);
x = x(:);

nData = length( x );

% Ground truth values
gt = {};
gt.k = 0.5;
gt.sigma = 1;

settings.dt = 0.1;
settings.nParticles = 1000;

dataOU = struct();
dataOU.x = data.X( 1, :);

ssFun = @(theta, data) -2 * ouParticleLikelihood( data, theta, settings);

model = struct();
model.ssfun  = ssFun;
model.sigma2 = 1.0;

model.N = nData;

nIterations = 1000;

optionsAM.nsimu = nIterations;
optionsAM.updatesigma = 0;
optionsAM.waitbar = 1;
optionsAM.verbosity = 0;
optionsAM.method = 'dram';

parameters = {
    {'Mean reversion rate: k', 0.5}
    {'Random walk noise level: sigma', 1, 0}
};

tic
[~, chainAM] = mcmcrun( model, dataOU, parameters, optionsAM);
toc

saveResults( "ou", chainAM, x);

