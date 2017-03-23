% This is a simple Interior Point Method to solve Convex QP for MPC
% Based on the Algorithm offered by the NTU Paper Embedded MPC and AUT's
% m-function "quad_wright".
% A simple version of priduip.(With centering path, without Pre.- Cor.)
% Yi Ding 2014.7.12

function [xstar,zstar,zstar1,Iter] = IPM_v2(G,GL,c,A,b,x,y,lambda)
global global_res sim_k;
[m,~] = size(A);
maxIteration = 30;
zstar1(1,1) = 0.5*x'*G*x + c'*x;
fx_old = 0.5*x'*G*x + c'*x;
x_old = x;
At = A';
sigma = 0.1;
Iter = 0;
%记录每次寻优的迭代周期值，测试一些参数的选择对寻优算法效率的影响；
Iter = 0;
global_res = zeros(maxIteration,1);
b_max = max(abs(b));
c_max = max(abs(c));

p_max = max(max(abs([G,c; A b])));

mu = y'*lambda/m;

%% 非迭代部分：检查解析解是否为全局最优
x_try = luEvaluate(GL,GL',-c);
if(min(A*x_try-b)>0)
    xstar = x_try;
    zstar = 0.5*x'*G*x + c'*x;
else
   
    for k = 1:maxIteration
        
        %% Section 1 Parameters Update
        rp_y = -A*x+b;
        rd = G*x-At*lambda+c;
        rp = -rp_y-y;
        
        inv_y_lambda = lambda./y;
        y_lambda_A = bsxfun(@times,A,inv_y_lambda);
        
        F = G+At*y_lambda_A;
         
        sig_mu = sigma*mu*ones(m,1);
        center_part = rp_y + sig_mu./lambda;
        equation_b = -rd+At*(inv_y_lambda.*center_part);
       
%        cond_F = cond(F)
        
        %% Section 2 Line Solve
        delta_x = line_solve(F,equation_b);
        delta_y = A*delta_x+rp;
        delta_lambda = -lambda-(lambda.*delta_y-sigma*mu*ones(m,1))./y;
        
        %% Section 3 Alpha Computing
        % AUT method to calculate alpha
        [alpha] = Alpha_AUT(y,lambda,delta_y,delta_lambda);

%         % NTU method to calculate alpha
%         [alpha] = Alpha_NTU(y,lambda,delta_y,delta_lambda);
        
        %% Section 4 Parameter Update
        mu_old = mu;
        
        x = x + alpha * delta_x;
        y = y + alpha * delta_y;
        lambda = lambda + alpha * delta_lambda;
        
        % Update Iterations
        Iter = Iter +1;
        
        mu = y'*lambda/m;
        %Solve centering parameter (Mehrotra's Heuristic)
        sigma = min((mu/mu_old)^3,0.99999);
        
        %% Section 5 Termination Check
%         % KKT residual (Wright) Termination Criteria
%         flag = TC_KKT_Wright(rd,rp,b,c,mu,0.00001);
%         if(flag==1)
%             break;        
%         end
        
%         % KKT residual (AUT) Termination Criteria
%         flag = TC_KKT_AUT(rd,rp,p_max,mu,0.00001);
%         if(flag==1)
%             break;        
%         end

%         % Duality Gap only (NTU)
%         if(mu<0.00001)
%             break;
%         end
% 
        % Himmelblau Termination Criteria
        [flag,fx_old,x_old] = TC_Himmelblau(G,c,A,b,fx_old,x_old,x,0.000001);
        if(flag==1)
            break;
        end
        


        
    end
xstar = x;
zstar = 0.5*x'*G*x + c'*x;
clc

end