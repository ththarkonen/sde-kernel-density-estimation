function [ x_t ] = ouPropagation( x0, theta, settings)
% this function is model specific
% Markov transition of the model f(x_{t} | x_{t-1})

    dt = settings.dt;
    nX = settings.nParticles;

    k = theta(1);
    sigma = theta(2);
    mu = 0;

    w_t = randn( nX, 1);
    x_t = x0 + k * ( mu - x0 ) * dt + sigma * sqrt(dt) * w_t;
end