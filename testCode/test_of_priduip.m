%This is a test file of priduip.m and a compare with Matlab function quadprog()
% X = priduip(G,c,A,b,x0,y0,lambda0) attempts to solve the quadratic programming 
%    problem:
%
%             min 0.5*x'*G*x + c'*x   subject to:  A*x >= b 
% Here x0,y0,lambda0 is the feasible start point.y0,lambda0 > 0;
% DingYi, 2013.10.17
% 
% H = [1,0;0,1];
% f = [1;1];
% A = [1,1;-1,2];
% b = [0;1]

%plot(0.5*x'*H*x + f'*x);

% [x,y]=meshgrid(linspace(0,2));
% %z = 0.5*H(1,1)*x^2 + 0.5*H(2,1)*x*y + 0.5*H(1,2)*x*y + 0.5*H(2,2)*y^2 + f(1,1)*x + f(2,1)*y
% z = x^2 + y^2;
% mesh(x,y,z);

clc;clear;
%According to the Poject in OR, the parameters are set as follows
G = [159.4797,144.1062;144.1062,134.9173];
c = [47.0497;43.0252];
A = [1,0;0,-1];
b = [0.4;0.4];

X = quadprog(G,c,A,b)
[star,zstar,zstar1] = priduip(G,c,-A,-b,[0;0],[1;1],[1;1])
[star,zstar,zstar1] = priduip_v2(G,c,-A,-b,[0;0],[1;1],[1;1])
%XXX = primal(G,c,A,b,[-2;-0.2],[0.001;0.1],[0.001;0.001])
%[xstar,zstar,zstar1] = barrier(G,c,A,b,[0;0]);

%stairs(zstar1,'linewidth',3)
%xlabel('Primal Dual Iteration');
%ylabel('f(x) value');
%title('f(x) value decrement');

