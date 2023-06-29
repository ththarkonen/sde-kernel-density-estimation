function [ll] = mertonParticleLikelihood( data, theta, settings)

    y = data.y;
    nData = length( y );

    ll = zeros( nData, 1);

    parfor ii = 2:nData

        y0 = y(ii - 1);
        y_ii = y(ii);

        [yPrediction] = mertonPropagation( y0, theta, settings);
        [yJumpPrediction] = mertonJumpPropagation( y0, theta, settings);

        p_ii = ksdensity( yPrediction, y_ii);
        pJump_ii = ksdensity( yJumpPrediction, y_ii);

        if( pJump_ii > p_ii )
            p_ii = pJump_ii;
        end

        ll_ii = log( p_ii );
        ll(ii) = ll_ii;
    end

    ll = sum( ll );
end