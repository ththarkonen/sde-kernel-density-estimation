function [ll] = ouParticleLikelihood( data, theta, settings)

    x = data.x;
    nData = length( x );

    ll = zeros( nData, 1);

    parfor ii = 2:nData

        x0 = x(ii - 1);
        x_ii = x(ii);

        xPrediction = ouPropagation( x0, theta, settings);
        density_ii = fitdist( xPrediction, 'Kernel');
        
        p_ii = pdf( density_ii, x_ii);
        ll_ii = log( p_ii );

        ll(ii) = ll_ii;
    end

    ll = sum( ll );
end