% Test the Kalman filter to prepare for real use
% 2014.8.22
% Yi

clc;clear;

load model2_wrong

nu = 1;ny = 2;nx = 4;   
N_sim  = 150;  

u = u_draw;
y = y_draw;
x_draw = x_draw';
xm_draw = xm_draw';

% Model2
Ac = [1 0.05  -0.0057  -0.000094883;
      0 1     -0.2308  -0.0057;
      0 0     1.0593   0.0510;
      0 0     2.3968   1.0593];
Bc = [0.0028;0.1106;-0.0091;-0.3674];
Cc = [1 0 0 0;
      0 0 1 0];
Dc = 0;

K_cons = [0.4624,-0.0189;0.4175,-0.3376;-0.0056,0.6674;-0.0489,3.8177;1.2928,-0.0209;-0.0076,1.5291];

[A_e,B_e,C_e] = augment(Ac,Bc,Cc,Dc);

R1 = eye(nx+ny,nx+ny);
R2 = eye(ny,ny);
P0 = eye(nx+ny,nx+ny);

x_est = zeros(N_sim,nx+ny);

% Here the noise in state and output is assumed to be independent
% Initial state x is assumed yo be zero here.
K{1} = (A_e*P0*C_e')*inv(C_e*P0*C_e'+R2);
x_est(1,:) = B_e*u(1)+K{1}*y(1,:)';
Sigma{1} = (A_e-K{1}*C_e)*P0*(A_e-K{1}*C_e)'+R1+K{1}*R2*K{1}';

% Riccati solve
[X_R,L_R,G_R] = dare(A_e',C_e',R1);
K_R = X_R*C_e'*inv(C_e*X_R*C_e'+R2);

for i=2:N_sim
   
    K{i} = (A_e*Sigma{i-1}*C_e')*inv(C_e*Sigma{i-1}*C_e'+R2);
    
    Sigma{i} = (A_e-K{i-1}*C_e)*Sigma{i-1}*(A_e-K{i}*C_e)'+R1+K{i}*R2*K{i}'; 
    
    %x_est(i,:) = A_e*x_est(i-1,:)'+B_e*u(i,:)'+K{i}*(y(i,:)'-C_e*x_est(i-1,:)');
    
    

    x_est(i,:) = A_e*x_est(i-1,:)'+B_e*u(i,:)'+K_cons*(y(i,:)'-C_e*x_est(i-1,:)');
    
end

x_draw_i = zeros(N_sim,4);
x_est_i = zeros(N_sim,4);
for i=2:N_sim
    x_draw_i(i,1:4) = x_draw(i,1:4) + x_draw(i-1,1:4);
    x_est_i(i,1:4) = x_est(i,1:4) + x_est(i-1,1:4);
end




