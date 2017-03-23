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
% V2 date: 2014.3.27 Modified after discussion with Kexin WANG
% V2 modification: linsolve scale: 3x3 to 1x1, delta_y,delta_lambda eliminate.


function [xstar,zstar,zstar1,Iter] = priduip_v2(G,c,A,b,x,y,lambda)
[m,n] = size(A);
maxIteration = 40;
zstar1(1,1) = 0.5*x'*G*x + c'*x;
res = 0.3;
tau = 0.7;
%记录每次寻优的迭代周期值，测试一些参数的选择对寻优算法效率的影响；
Iter = 0;

for k = 1:maxIteration
    
    rd = G*x-A'*lambda+c;
    rp = A*x-y-b;
    
    diag_y = diag(y);
    inv_diag_y = diag_inv(diag_y);
    diag_lambda = diag(lambda);
    inv_diag_lambda = diag_inv(diag_lambda);

    
    F = G+A'*inv_diag_y*diag_lambda*A;
    equation_b = -rd+A'*inv_diag_y*diag_lambda*(-rp-y);
    %Left = linsolve(F,equation_b);
    Left = line_solve(F,equation_b);
    delta_x_aff = zeros(n,1);
    delta_y_aff = zeros(m,1);
    delta_lambda_aff = zeros(m,1);
    for kk = 1:n
        delta_x_aff(kk,:) = Left(kk,:);
    end
    delta_y_aff = A*delta_x_aff+rp;
    delta_lambda_aff = -lambda-inv_diag_y*diag_lambda*delta_y_aff;
    
    mu = y'*lambda/m; 
    
    alpha_aff = 1;
    while(max(((y+alpha_aff*delta_y_aff)<0)|((lambda+alpha_aff*delta_lambda_aff)<0)))   %Attention!&&~
        alpha_aff = alpha_aff-0.01;        %每次减半？
        if(alpha_aff<=0)
            break;
        end
    end
        
    mu_aff = (y + alpha_aff * delta_y_aff)'*(lambda + alpha_aff * delta_lambda_aff)/m;
    
    sigma = (mu_aff/mu)^3;
    
    equation_b = -rd+A'*inv_diag_y*diag_lambda*(-rp-y+sigma*mu*inv_diag_lambda*ones(m,1)-inv_diag_lambda*diag(delta_lambda_aff)*diag(delta_y_aff)*ones(m,1));
    delta_x = zeros(n,1);
    delta_y = zeros(m,1);
    delta_lambda = zeros(m,1);
    %Left = linsolve(F,equation_b);
    Left = line_solve(F,equation_b);
    for kk = 1:n
        delta_x(kk,:) = Left(kk,:);
    end
    delta_y = A*delta_x+rp;
    delta_lambda = -lambda-inv_diag_y*(diag_lambda*delta_y+diag(delta_lambda_aff)*diag(delta_y_aff)*ones(m,1)-sigma*mu*ones(m,1));
    
    %迭代快结束时，使tau值趋近于1；
    if(res<0.3)
        tau = 1-res;
    end
    
    alpha_tau_pri = 1;
    temp_cond = ((y+alpha_tau_pri*delta_y)<(1-tau)*y);
    while(sum(temp_cond)>0)
        alpha_tau_pri = alpha_tau_pri-0.01;    %每次减半？
        if(alpha_tau_pri<=0)
            break;
        end
        temp_cond = ((y+alpha_tau_pri*delta_y)<(1-tau)*y);
    end
    alpha_tau_dual = 1;
    temp_cond = ((lambda+alpha_tau_dual*delta_lambda)<(1-tau)*lambda);
    while(sum(temp_cond)>0)
        alpha_tau_dual = alpha_tau_dual-0.01;   %每次减半？
        if(alpha_tau_dual<=0)
            break;
        end
        temp_cond = ((lambda+alpha_tau_dual*delta_lambda)<(1-tau)*lambda);
    end
    alpha = min(alpha_tau_pri,alpha_tau_dual);
    
    x = x + alpha * delta_x;
    y = y + alpha * delta_y;
    lambda = lambda + alpha * delta_lambda;
    
    %计算KKT残差
    res_1 = abs(G*x-A'*lambda+c);
    temp1 = G*x;
    temp2 = A'*lambda-c;
    temp2 = temp1-temp2;
    res_2 = abs(A*x-y-b);
    res_3 = abs(mu);
    res = max([res_1;res_2;res_3]);
    
    %迭代次数加1
    Iter = Iter +1;
%     disp('The Iteration is:');
%     disp(Iter);
%     disp('The Residual is:');
%     disp(res);
    
%     zstar1(k+1,1) = 0.5*x'*G*x + c'*x;
    
    %终止条件
    if(res<0.001)
        break;
    end        
end
xstar = x;
zstar = 0.5*x'*G*x + c'*x;
end