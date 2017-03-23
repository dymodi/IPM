%反馈矫正环节
%这里对全量状态进行预测和校正
%2014.8.22
%还是不大对，还是直接观测增量状态吧

function [x_k] = fankui_v3(x_k_old,y_k,u_k_1,u_k_2)
global A_e B_e C_e;

global R1 R2 P0;
global K x_est Sigma;


% Kalman matrix (Created by the routine kal_filter).
% This matrix is only designed for model2.
K = [0.4624,-0.0189;0.4175,-0.3376;-0.0056,0.6674;-0.0489,3.8177;1.2928,-0.0209;-0.0076,1.5291];

% %预测
% x_k = A_e * x_k_1 + B_e * (u_k_1 - u_k_2); 
% %矫正
% x_k = x_k + L * (y_k - C_e * x_k); 


    K = (A_e*Sigma*C_e')*inv(C_e*Sigma*C_e'+R2);
    Sigma = (A_e-K*C_e)*Sigma*(A_e-K*C_e)'+R1+K*R2*K'; 
    
    x_est = A_e*x_est+B_e*(u_k_1-u_k_2)+K*(y_k-C_e*x_est);
    
    
    
    x_k = x_est;


end