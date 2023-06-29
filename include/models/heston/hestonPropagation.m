function [ SV_t] = hestonPropagation( S0, V0, parameters, settings)
% this function is model specific
% Markov transition of the model f(x_{t} | x_{t-1})

    r = settings.r;
    q = settings.q;

    T = settings.dt;
    nSteps = settings.simulationSteps;
    nX = settings.nParticles;
    
    dT = T / nSteps;

    theta = parameters(1);
    kappa = parameters(2);
    xi = parameters(3);
    rho = parameters(4);
    lambda_1 = parameters(5);
    
    eta = lambda_1 * sqrt( 1 - rho^2);
    
    lnS = zeros( nX, nSteps + 1);
    V = zeros( nX, nSteps + 1);
    
    lnS( :, 1) = log( S0 );
    V( :, 1) = V0;
    
    b = eta -.5;
    
    k1 = exp(-kappa*dT);
    k2 = xi^2*k1.*(1-k1)/kappa;
    k3 = exp(kappa*dT)*0.5.*k2.*(1-k1).*theta;

    phiC = 1.5;                             
    gamma1 = .5;                            
    gamma2 = .5;                            

    c1 = (r-q)*dT;                          
    c2 = -rho*kappa/xi*dT;             

    K0 = c1 + c2.*theta;                    
    K1 = gamma1*dT*(kappa*rho/xi+b)-rho/xi; 
    K2 = gamma2*dT*(kappa*rho/xi+b)+rho/xi; 
    K3 = gamma1*dT*(1-rho^2);                            
    K4 = gamma2*dT*(1-rho^2);                            

    UV1 = rand( nX, nSteps);         
    UV = rand( nX, nSteps);         
    dW2 = norminv( UV );             

    for ii = 2:( nSteps + 1 )
        
        m = theta + (V(:,ii-1)-theta)*k1; 
        s2 = V(:,ii-1)*k2 + k3;           
        phi = s2./m.^2;                                  
        phihat = 1./phi;                                    
        b2 = 2*phihat - 1 + sqrt(2*phihat.*(2*phihat-1));   
        a = m ./ (1 + b2);                                      
        
        if isempty(V(phi<=phiC,ii))
        else
            V(phi<=phiC,ii) = a(phi<=phiC).*(sqrt(b2(phi<=phiC)) + norminv(UV1(phi<=phiC,ii-1))).^2;                
        end
        
        p = (phi - 1)./(phi + 1);
        V((UV1(:,ii-1)<=p) & (phi>phiC),ii) = 0; 
        I1b = find((UV1(:,ii-1)>p) & (phi>phiC));    
        beta = (1 - p)./m;  
        
        if isempty(I1b)
        else       
            V(I1b,ii) = log((1-p(I1b))./(1-UV1(I1b,ii-1)))./beta(I1b);  
        end
        
        lnS(:,ii) = lnS(:,ii-1) + K0 + K1.*V(:,ii-1) + K2.*V(:,ii) ... 
                + sqrt(K3.*V(:,ii-1) + K4.*V(:,ii)).*dW2(:,ii-1);
    end
    
    lnS_t = lnS( :, end);
    S_t = exp( lnS_t ); 
    V_t = V( :, end);
    
%     pathS = exp( lnS );
%     pathV = V;
    SV_t = [ S_t, V_t];
end