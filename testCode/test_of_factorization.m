% Test of factorization 
% This is a m file to test all kinds of matrix factorization
% DingYi
% 2013.12.5

clc;clear;
% Test of Cholesky Factorization
G = [43.9736,32.0645;
     32.0645,26.1758];
A = [-1,0;0,-1;1,0;0,1;-0.5,0;-1.75,-0.5];
lambda_y = diag([10,10,3,3,7.5,6]);
equation_b = [-22.5755;-20.2424];

B = A'*lambda_y*A;
F = G+B;

F = [4.6460,3.8305,3.3318;
     3.8305,3.7666,3.0105;
     3.3318,3.0105,3.0020];

% L1 = cf(G)
% L2 = cf(B)
% A = [0.4694    0.1622    0.5285
%      0.0119    0.7943    0.1656
%      0.3371    0.3112    0.6020]
% tic
% for k = 1:1000
L1 = cf(F)
% end
L1*L1'
% toc
% 
% tic
% for k = 1:1000
% [L2,D2] = cf_v2(F);
% end
% L2*D2*L2'
% toc

% tic
% for k = 1:1000
 [L3,D3] = mcf(F)
% end
L3*D3*L3'
L3*D3
% toc
% Left = linsolve(F,equation_b);
% luEvaluate(L,L',equation_b);

% L = [8.7877,0;4.2462,4.7587];
% b = 

% Test of Symmetric Indefinite Factorization
% A1 =[1     2     3
%      2     1     3
%      3     3     2];
% A2 =[1,2,3,4;2,3,4,5;3,4,5,6;4,5,6,7];
% K =[
%      6     2     1     1     0
%      2     5     2     0     1
%      1     2     4     1     1
%      1     0     1     0     0
%      0     1     1     0     0];
% % [P,L,D] = sif1(A1)
% % L*D*L'
% [L,D] = ldl(A2)
% [P,L,D] = sif2(A2)
% P*A2*P'
% L*D*L'
