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
% V3 date: 2014.5.20 Modified after discussion with Zhao and Xu. Confirm
% the priduip as function.
% V3 modification: y and lambda computing,alpha decrease,
% V4 date: 2014.5.27 Modified after discussion with Shao, Xu and Chen. 
% V4 modification: residuel meaning.



function [xstar,zstar,zstar1,Iter] = priduip_v4(G,GL,c,A,b,x,y,lambda)
global global_res;
[m,~] = size(A);
maxIteration = 30;
zstar1(1,1) = 0.5*x'*G*x + c'*x;
fx_old = 0.5*x'*G*x + c'*x;
x_old = x;
zstar = 0;
res = 0.3;
tau = 0.7;
At = A';
%记录每次寻优的迭代周期值，测试一些参数的选择对寻优算法效率的影响；
Iter = 0;
global_res = zeros(maxIteration,1);
b_max = max(abs(b));
c_max = max(abs(c));

p_max = max(max(abs([G,c; A b])));

mu = y'*lambda/m;

% %检查解析解是否为全局最优
x_try = luEvaluate(GL,GL',-c);
if(min(A*x_try-b)>=0)
    xstar = x_try;
else
for k = 1:maxIteration
    
    rp_y = -A*x+b;
    rd = G*x-At*lambda+c;
    rp = -rp_y-y;
    
    inv_y_lambda = lambda./y;
    y_lambda_A = bsxfun(@times,A,inv_y_lambda);
    
    F = G+At*y_lambda_A;
    equation_b = -rd+At*(inv_y_lambda.*rp_y);   
    delta_x_aff = line_solve(F,equation_b);
    delta_y_aff = A*delta_x_aff+rp;
    delta_lambda_aff = -lambda-inv_y_lambda.*delta_y_aff;
    
    mu = y'*lambda/m; 
    
    %Decide on suitable affine step-length
    duals = [y;lambda];
    delta = [delta_y_aff;delta_lambda_aff];
    index = delta < 0; 
    if(any(index))
        alpha_aff = 0.9995*min(-duals(index)./delta(index)); %solves for min ratio (max value of alpha allowed)
    else
        alpha_aff = 0.999995;
    end
    
    %Check for numerical problems (alpha will go negative)
    if(alpha_aff < 0), error('alpha<0'); end   
    
    mu_aff = (y + alpha_aff * delta_y_aff)'*(lambda + alpha_aff * delta_lambda_aff)/m;
    
    sigma = (mu_aff/mu)^3;
    
    sig_mu_lam_y = sigma*mu*ones(m,1)-delta_lambda_aff.*delta_y_aff;
    aff_part = rp_y + sig_mu_lam_y./lambda;
    equation_b = -rd+At*(inv_y_lambda.*aff_part);
    %equation_b = -rd+A'*inv_diag_y*diag_lambda*(-rp-y+sigma*mu*inv_diag_lambda*ones(m,1)-inv_diag_lambda*diag(delta_lambda_aff)*diag(delta_y_aff)*ones(m,1));
    
    %Left = linsolve(F,equation_b);
    delta_x = line_solve(F,equation_b);
    delta_y = A*delta_x+rp;
    %delta_lambda = -lambda-inv_diag_y*(diag_lambda*delta_y+diag(delta_lambda_aff)*diag(delta_y_aff)*ones(m,1)-sigma*mu*ones(m,1));
    delta_lambda = -lambda-(lambda.*delta_y+delta_lambda_aff.*delta_y_aff-sigma*mu*ones(m,1))./y;
    
    %迭代快结束时，使tau值趋近于1；
    if(res<0.3)
        tau = 1-res;
    end
    
    index = delta_y < 0; 
    if(any(index))
        alpha_tau_pri = 0.9995*min(-tau*y(index)./delta_y(index)); %solves for min ratio (max value of alpha allowed)
    else
        alpha_tau_pri = 0.999995;
    end
    %Check for numerical problems (alpha will go negative)
    if(alpha_tau_pri < 0), error('alpha<0'); end  
            
    index = delta_lambda < 0; 
    if(any(index))
        alpha_tau_dual = 0.9995*min(-tau*lambda(index)./delta_lambda(index)); %solves for min ratio (max value of alpha allowed)
    else
        alpha_tau_dual = 0.999995;
    end
    %Check for numerical problems (alpha will go negative)
    if(alpha_tau_dual < 0), error('alpha<0'); end  
     
    alpha = min(alpha_tau_pri,alpha_tau_dual);
    
    x = x + alpha * delta_x;
    y = y + alpha * delta_y;
    lambda = lambda + alpha * delta_lambda;
  
    % Update Iterations
    Iter = Iter +1;
%   disp('The Iteration is:');
%   disp(Iter);
%   disp('The Residual is:');
%   disp(res);
    
%     % KKT residuel and Duality measure Criteria 
%     res_1 = abs(G*x-A'*lambda+c);
%     res_2 = abs(A*x-y-b);
%     res_3 = abs(mu);
%     res = max([res_1;res_2;res_3]);
%     %终止条件
%     if(res<0.00001)
%         break;
%     end 

%     % KKT residual (AUT) Termination Criteria
%     flag = TC_KKT_AUT(rd,rp,p_max,mu,0.0001);
%     if(flag==1)
%         break;
%     end

%     % KKT residual (Wright) Termination Criteria
%     flag = TC_KKT_Wright(rd,rp,b,c,mu,0.0001);
%     if(flag==1)
%         break;
%     end

%     % KKT residuel and Duality measure Criteria based on PD-IPM p226
%     res_1 = max(abs(G*x-A'*lambda+c))/(1+c_max);
%     res_2 = max(abs(A*x-y-b))/(1+b_max);
%     res_3 = abs(c'*x-b'*y)/(1+abs(c'*x));
%     res = max([res_1;res_2;res_3]);
%     % Check residuel
%     if(res<0.00001)
%         break;
%     end   
    
    % Himmelblau Termination Criteria
    fx = 0.5*x'*G*x + c'*x; 
    res1 = abs(fx_old - fx);
    res2 = abs(x_old - x);
    cons_obey = min(A*x-b);
    res = max([res1;res2]);
    global_res(Iter+1) = res;
    fx_old = fx;
    x_old = x;
    % Check residuel
    if(cons_obey>=-0.0001&&res<0.0001)
        break;
    end 
    
    

end
xstar = x;
zstar = 0.5*x'*G*x + c'*x;

% %Draw for RPS preparation.
% if(Iter>9)
%     figure;
%     plot([Iter-4:Iter],global_res(Iter-4:Iter),'linewidth',2);
%     xlabel('QP Iteration');
%     ylabel('Residuel');
% end

% end

end