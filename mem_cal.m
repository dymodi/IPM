% Calculate the estimation of the memory

clc;clear;

% %Model1
% M = 10;P = 20;
% nx = 2;ny = 1;nyc = 0;nu = 1;
% ndec = nu*M;
% mc = 4*nu*M + 2*nyc*P;

% % Model2
% M = 10;P = 25;
% nx = 4;ny = 2;nyc = 2;nu = 1;
% ndec = nu*M;
% mc = 2*nu*M + 2*nyc*P;

%Model3
M = 10;P = 10;
nx = 4;ny = 3;nyc = 1;nu = 1;
ndec = nu*M;
mc = 4*nu*M + 2*nyc*P;

%// Memory allocation for main function
M_main = (nx+ny)*ny+P*ny*(nx+ny)+M*nu*P*ny+P*nu*M*nu*M+mc*ndec+mc*ndec+ndec*ndec*3+ndec+mc*2+P*ny;

%// Memory allocation for mpc_DP
M_mpc = ndec*4+mc+ny*P+ny*P*M*nu+max([ny*P*(nx+ny),mc,ndec*ndec,mc*ndec]);

%//Memory allocation for IPM_V3_DP
M_IPM = ndec*5+mc*7+mc*ndec+ndec*ndec+mc*2+max([mc,ndec*ndec]);

%//Memory allocation for feedback
M_feed = 2*(nx+ny);

%//Memory allocation for mcf
M_mcf = ndec*ndec+ndec+2;

%//Memory allocation for alpha_decreas
M_alpha = mc*4;

%//Memory allocation for line_solve
M_line = ndec*ndec*4+ndec*3;

M_total = M_main + M_mpc + M_IPM + M_feed + M_mcf + M_alpha + M_line;

KB_total = M_total*8/1024;
