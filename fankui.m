%反馈矫正环节
function [x_k] = fankui(flag,x_k_1,y_k,u_k_1,u_k_2,A_e,B_e,C_e)
L=[0;0;1];  

if(flag==1)
    x_k = [0;0;y_k];
else
    %预测
    x_k = A_e * x_k_1 + B_e * (u_k_1 - u_k_2); 
    %矫正
    x_k = x_k + L * (y_k - C_e * x_k); 
end