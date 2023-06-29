function [ x_t ] = mertonPropagation( x0, theta, parameters)
% this function is model specific
% Markov transition of the model f(x_{t} | x_{t-1})

    dt = parameters.dt;
    nX = parameters.nParticles;

    k = theta(1);
    sigma = theta(2);
    jumpMu = theta(3);
    jumpDelta = theta(4);

    w_t = randn( nX, 1);
    x_t = x0 - k * x0 * dt + sigma * sqrt(dt) * w_t;
end