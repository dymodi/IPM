% This is a QP solving function based on Primal-dual Interior Point Method
% Based on Algorithm 16.4 Predictor-Corrector Algorothm for QP in book‘Numerical Optimization’
% X = priduip(G,c,A,b,x0,y0,lambda0) attempts to solve the quadratic programming 
%    problem:
%             min 0.5*x'*G*x + c'*x   subject to:  A*x >= b 
% Here x0,y0,lambda0 is the feasible start point.y0,lambda0 > 0;
% DingYi,
% Start: 2013.10.17
% Debug: 2013.10.18 a whole day, still not correct.
% Successful tested: 2013.10.19

function [xstar,zstar,zstar1] = priduip(G,c,A,b,x,y,lambda)
[m,n] = size(A);
maxIteration = 40;
tau = zeros(maxIteration,1);
zstar1(1,1) = 0.5*x'*G*x + c'*x;

for k = 1:maxIteration
    
    rd = G*x-A'*lambda+c;
    rp = A*x-y-b;
    
    F = cell2mat({G,zeros(n,m),-A';A,-eye(m,m),zeros(m,m);zeros(m,n),diag(lambda),diag(y);});
    equation_b = [-rd;-rp;-diag(lambda)*diag(y)*ones(m,1)];
    Left = linsolve(F,equation_b);
    delta_x_aff = zeros(n,1);
    delta_y_aff = zeros(m,1);
    delta_lambda_aff = zeros(m,1);
    for kk = 1:n
        delta_x_aff(kk,:) = Left(kk,:);
    end
    for kk = n+1:n+m
        delta_y_aff(kk-n,:) = Left(kk,:);
    end
    for kk = n+m+1:n+2*m
        delta_lambda_aff(kk-n-m,:) = Left(kk,:);
    end
    
    mu = y'*lambda/m; 
    
    alpha_aff = 1;
    while(max(((y+alpha_aff*delta_y_aff)<0)|((lambda+alpha_aff*delta_lambda_aff)<0)))   %Attention!&&~
        alpha_aff = alpha_aff - 0.01;
        if(alpha_aff<=0)
            break;
        end
    end
        
    mu_aff = (y + alpha_aff * delta_y_aff)'*(lambda + alpha_aff * delta_lambda_aff)/m;
    
    sigma = (mu_aff/mu)^3;
    
    equation_b = [-rd;-rp;-diag(lambda)*diag(y)*ones(m,1)-diag(delta_lambda_aff)*diag(delta_y_aff)*ones(m,1)+sigma*mu*ones(m,1)];
    delta_x = zeros(n,1);
    delta_y = zeros(m,1);
    delta_lambda = zeros(m,1);
    Left = linsolve(F,equation_b);
    for kk = 1:n
        delta_x(kk,:) = Left(kk,:);
    end
    for kk = n+1:n+m
        delta_y(kk-n,:) = Left(kk,:);
    end
    for kk = n+m+1:n+2*m
        delta_lambda(kk-n-m,:) = Left(kk,:);
    end

    tau(k,1) = 0.6;      %Atten! tau(k) should be set to approach to 1 when iterates approach solution.
    alpha_tau_pri = 1;
    temp_cond = ((y+alpha_tau_pri*delta_y)<(1-tau(k,1))*y);
    while(sum(temp_cond)>0)
        alpha_tau_pri = alpha_tau_pri - 0.01;    %每次减半
        if(alpha_tau_pri<=0)
            break;
        end
        temp_cond = ((y+alpha_tau_pri*delta_y)<(1-tau(k,1))*y);
    end
    alpha_tau_dual = 1;
    temp_cond = ((lambda+alpha_tau_dual*delta_lambda)<(1-tau(k))*lambda);
    while(sum(temp_cond)>0)
        alpha_tau_dual = alpha_tau_dual - 0.01;   %每次减半
        if(alpha_tau_dual<=0)
            break;
        end
        temp_cond = ((lambda+alpha_tau_dual*delta_lambda)<(1-tau(k))*lambda);
    end
    alpha = min(alpha_tau_pri,alpha_tau_dual);
    
    x = x + alpha * delta_x;
    y = y + alpha * delta_y;
    lambda = lambda + alpha * delta_lambda;
    zstar1(k+1,1) = 0.5*x'*G*x + c'*x;
    
    %if(mu<0.00001)                           %end criterion for mu
    if(abs(zstar1(k+1,1)-zstar1(k,1))<1e-8)    %end criterion for TolFun，
        break;
    end        
end
xstar = x;
zstar = 0.5*x'*G*x + c'*x;
end