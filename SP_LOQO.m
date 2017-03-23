% This is one of the ways to choose a starting point -- LOQO's method.
% Based on the file LOQO: AN INTERIOR POINT CODE FOR QUADRATIC PROGRAMMING
% 2014.7.13
% Yi

function [delta_u_ini,y_ini,lambda_ini] = SP_LOQO(G,c,A,b)

At = A';

[m,n] = size(G);
F = G+eye(m,n)+At*A;
equation_b = At*b-c;
delta_u_ini = line_solve(F,equation_b);
lambda_ini = b-A*delta_u_ini;

[m,~] = size(b);
y_ini = max(abs(lambda_ini),1*ones(m,1));


end