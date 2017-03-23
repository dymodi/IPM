function [delta_u,u_k,y_k,x_k] = mpcTime(delta_u,u_k,y_k,x_k,r_k)

global Ac Bc Cc nu ny P M;
global Phi F Q xm G;
global u_k_1 u_k_2;
global delta_U_p delta_U_n U_p U_n Y_p Y_n OMEGA_L;

x_k = fankui_v2(x_k,y_k,u_k_1,u_k_2);
aug_u_k_1 = [];
for k = 1:M
    aug_u_k_1 = [aug_u_k_1;u_k_1];
end
%带被控变量的约束
omega_r = [delta_U_p;-delta_U_n;U_p-aug_u_k_1;-U_n+aug_u_k_1;Y_p-F*x_k;-Y_n+F*x_k];
%不带被控变量的约束
%omega_r = [delta_U_p;-delta_U_n;U_p-aug_u_k_1;-U_n+aug_u_k_1];
c = (F*x_k-r_k)'*(Q*eye(ny*P,ny*P))*Phi;
c = c';
%delta_u = quadprog(G,c,OMEGA_L,omega_r,[],[],[],[],[0;0],'TolFun',1e-8);

%确定每个周期优化求解的初值
delta_u_ini = [];
for k = 1:M
    delta_u_ini = [delta_u_ini;delta_u];
end
[m,~] = size(omega_r);
y_ini = zeros(m,1);
h_x = -OMEGA_L*delta_u_ini+omega_r;
for i=1:m
y_ini(i) = max([h_x(i),0.5]);
end;
lambda_ini = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;0.1;0.2;0.3;0.4;0.5;0.6];

%内点法求解优化命题
[delta_u,~,~] = priduip(G,c,-OMEGA_L,-omega_r,delta_u_ini,y_ini,lambda_ini);
%[delta_u,~,~,Iter] = priduip_v2(G,c,-OMEGA_L,-omega_r,delta_u_ini,y_ini,lambda_ini);
delta_u = delta_u(1:2,1);
u_k = u_k + delta_u;
end