
function [delta_u,delta_u_M_out,u_k,y_k,x_k,delta_u_ini,y_ini,lambda_ini] = mpc_v3_model2(delta_u,delta_u_M_in,u_k,y_k,x_k,r_k,delta_u_ini,y_ini,lambda_ini)

global Ac Bc Cc nu ny P M sim_k A_e B_e nx;
global Phi Phi1 Phi2 F F1 F2 Q QQ QQ2 xm G GL xr xm_old;
global u_k_1 u_k_2;
global delta_U_p delta_U_n U_p U_n Y_p Y_n OMEGA_L;
global QPTime TimeIter TotalIter;

global yr yp;

global xm2;

% % Feedback
% x_k = fankui_v2(x_k,y_k,u_k_1,u_k_2);

%V3 uses Kalman filter
x_k = fankui_v3(x_k,y_k,u_k_1,u_k_2);
% 
% % Direct use
% x_k = xr;

aug_u_k_1 = [];
for k = 1:M
    aug_u_k_1 = [aug_u_k_1;u_k_1];
end
% %带被控变量的约束
% omega_r = [delta_U_p;-delta_U_n;U_p-aug_u_k_1;-U_n+aug_u_k_1;Y_p-F*x_k;-Y_n+F*x_k];
%带被控变量不带控制增量的约束
omega_r = [U_p-aug_u_k_1;-U_n+aug_u_k_1;Y_p-F*x_k;-Y_n+F*x_k];
%不带被控变量的约束
%omega_r = [delta_U_p;-delta_U_n;U_p-aug_u_k_1;-U_n+aug_u_k_1];
c = (F*x_k-r_k)'*QQ*Phi;
c = c';
%c = c*10^-35;

%确定初值的DY的方法
[delta_u_ini,y_ini,lambda_ini] = SP_DY(delta_u_M_in,omega_r,OMEGA_L,0.5,delta_u_ini,y_ini,lambda_ini);

%确定初值的AUT的方法
%[delta_u_ini,y_ini,lambda_ini] = SP_AUT(G,c,-OMEGA_L,-omega_r,0.5);
% 
% %确定初值的LOQO的方法（有问题，暂时搁置）
%[delta_u_ini,y_ini,lambda_ini] = SP_LOQO(G,c,-OMEGA_L,-omega_r);
% 
% %确定初值的Wright的启发式方法（需要用户提供一个初始点）
% [delta_u_ini,y_ini,lambda_ini] = SP_wright(delta_u_ini,y_ini,lambda_ini,G,c,-OMEGA_L,-omega_r);



tic
%内点法求解优化命题
%[delta_u_M_out,~,~,Iter] = priduip_v2(G,c,-OMEGA_L,-omega_r,delta_u_ini,y_ini,lambda_ini);
%[delta_u_M_out,~,~,Iter] = priduip_v4(G,GL,c,-OMEGA_L,-omega_r,delta_u_ini,y_ini,lambda_ini);
[delta_u_M_out,~,~,Iter] = IPM_v2(G,GL,c,-OMEGA_L,-omega_r,delta_u_ini,y_ini,lambda_ini);
%[delta_u,~,Iter,~,~] = quad_wright(G,c,OMEGA_L,omega_r,60,0.00001,0,delta_u_ini,lambda_ini,y_ini);
QPTime(TimeIter,1) = toc;
TotalIter(TimeIter,1) = Iter;
TimeIter = TimeIter + 1;
%简易内点法
%[delta_u,~,~,Iter,y,lambda] = ipm(G,c,-OMEGA_L,-omega_r,delta_u_ini,y_ini,lambda_ini);
%简易热启动内点法


delta_u = delta_u_M_out(1:nu,1);

xm_old = xm;

u_k = u_k + delta_u;

% %Linear Model
% xm = Ac*xm + Bc*u_k;  %全量状态更新
% y_k = Cc*xm;            %输出更新

%Nonlinear Model
[y_k,xm] = MPConDSP_mode2(u_k,xm);

u_k_2 = u_k_1;
u_k_1 = u_k;

% Acquire x_r;2014.8.22 Add
xr = A_e * xr + B_e * delta_u;
yp = F * xr + Phi * delta_u_M_out;
xr(1:nx) = xm - xm_old;
xr(nx+1:nx+ny) = y_k;

end