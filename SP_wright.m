% This is one of the ways to choose a starting point -- Wright's heuristic
% method.
% Based on the book Numerical Optimization, P484
% 2014.7.13
% Yi

function [delta_u_ini,y_ini,lambda_ini] = SP_wright(delta_u_bar,y_bar,lambda_bar,G,c,A,b)

At = A';

rp_y = -A*delta_u_bar+b;
rd = G*delta_u_bar-At*lambda_bar+c;
rp = -rp_y-y_bar;

inv_y_lambda = lambda_bar./y_bar;
y_lambda_A = bsxfun(@times,A,inv_y_lambda);

F = G+At*y_lambda_A;
equation_b = -rd+At*(inv_y_lambda.*rp_y);
delta_x_aff = line_solve(F,equation_b);
delta_y_aff = A*delta_x_aff+rp;
delta_lambda_aff = -lambda_bar-inv_y_lambda.*delta_y_aff;

y_ini = max(1,abs(y_bar+delta_y_aff));
lambda_ini = max(1,abs(lambda_bar+delta_lambda_aff));
delta_u_ini = delta_u_bar;

end
