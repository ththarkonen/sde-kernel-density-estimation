clc
clear
close all

includeFolders = genpath('include');
addpath( includeFolders );

filePaths = ["data/merton/merton1.mat",...
             "data/merton/merton2.mat",...
             "data/merton/merton3.mat"];

nFiles = length( filePaths );

for ii = 1:nFiles

    filePath_ii = filePaths(ii)

    data = load( filePath_ii );
    nData = length( data.y );
    
    settings.dt = 0.01;
    settings.nParticles = 15000;
    
    dataMerton = struct();
    dataMerton.y = data.y;
    
    ssFun = @(theta, data) -2 * mertonParticleLikelihood( data, theta, settings);
    
    model = struct();
    model.ssfun  = ssFun;
    model.sigma2 = 1.0;
    
    model.N = nData;
    nIterations = 50000;
    
    optionsDRAM.nsimu = nIterations;
    optionsDRAM.updatesigma = 0;
    optionsDRAM.waitbar = 0;
    optionsDRAM.verbosity = 2;
    optionsDRAM.method = 'dram';
    optionsDRAM.burnintime = 100;
    optionsDRAM.adaptint = 100;
    optionsDRAM.drscale = 5;
    
    parameters = {
        {'k', 10}
        {'sigma', 0.08, 0}
        {'jumpMu', 0.01}
        {'jumpDelta', 0.1, 0}
    };
    
    tic
    [~, chainAM] = mcmcrun( model, dataMerton, parameters, optionsDRAM);
    toc
    saveResults( "merton", chainAM, data.y);
end
