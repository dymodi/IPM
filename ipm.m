% This is a simple Interior Point Method to solve Convex QP for MPC
% Based on the Algorithm offered by the NTU Paper Embedded MPC
% A simple version of priduip.
% Yi Ding 2014.5.4

function [xstar,zstar,zstar1,Iter,y,lambda] = ipm(G,c,A,b,x,y,lambda)
[m,n] = size(A);
maxIteration = 10;
zstar1(1,1) = 0.5*x'*G*x + c'*x;
sigma = 0.25;
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
    %equation_b = -rd+A'*inv_diag_y*diag_lambda*(-rp-y);
    equation_b = -rd+A'*inv_diag_y*diag_lambda*(-rp-y+sigma*mu*inv_diag_lambda*ones(m,1));
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
    
    alpha_aff = 1;
    while(max(((y+alpha_aff*delta_y_aff)<0)|((lambda+alpha_aff*delta_lambda_aff)<0)))   %Attention!&&~
        alpha_aff = alpha_aff-0.01;        %每次减半？
        if(alpha_aff<=0)
            break;
        end
    end
    
    x = x + alpha_aff * delta_x_aff;
    y = y + alpha_aff * delta_y_aff;
    lambda = lambda + alpha_aff * delta_lambda_aff;
    
    mu = y'*lambda/m; 
    
    %计算KKT残差
    res_1 = abs(G*x-A'*lambda+c);
    res_2 = abs(A*x-y-b);
    res_3 = abs(mu);
    res = max([res_1;res_2;res_3]);
    
    %迭代次数加1
    Iter = Iter +1;
%     disp('The Iteration is:');
%     disp(Iter);
%     disp('The /mu is:');
%     disp(mu);
    
    %终止条件
    if(res<0.001)
        break;
    end    
end
xstar = x;
zstar = 0.5*x'*G*x + c'*x; 
end
