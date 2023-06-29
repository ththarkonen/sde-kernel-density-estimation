function [ll] = hestonParticleLikelihood( data, theta, settings)

    S = data.S;
    V = data.V;
    nData = length( S );

    ll = zeros( nData, 1);

    flag = constraints( theta );

    if( flag == 0 )

        ll = -Inf;
        return;
    end

    parfor ii = 2:nData

        S0 = S(ii - 1);
        S_ii = S(ii);
        
        V0 = V(ii - 1);
        V_ii = V(ii);
        
        SV_ii = [S_ii, V_ii];

        [svPrediction] = hestonPropagation( S0, V0, theta, settings);
        
        p_ii = ksdensity( svPrediction, SV_ii);
        ll_ii = log( p_ii );

        ll(ii) = ll_ii;
    end

    ll = sum( ll );
end

function flag = constraints( parameters )

    theta = parameters(1);
    kappa = parameters(2);
    xi = parameters(3);
    rho = parameters(4);
    lambda_1 = parameters(5);
    
    flag = 1;
    
    if any( [parameters(1:3)] < 0 )
        flag = 0;
    end  %positivity

    fellerCondition = 2 * theta * kappa - xi^2;

    % Feller condition
    if( fellerCondition < 0)
        flag = 0;
    end

    % correlation
    if abs( rho ) > 1
        flag = 0; 
    end
end
