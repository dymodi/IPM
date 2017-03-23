% This is one of the ways to choose a starting point -- DY's method.
% Based on the book Numerical Optimization, P484
% 2014.7.13
% Yi

function [delta_u_ini,y_ini,lambda_ini] = SP_DY(delta_u_M_in,omega_r,OMEGA_L,WarmVal,delta_u_ini,y_ini,lambda_ini)

global nu M;

%delta_u starting point 
delta_u_ini = zeros(nu*M,1);
for k=1:nu*(M-1)
    delta_u_ini(k) = delta_u_M_in(k+nu);
end

% % Easy way
% lambda_ini = WarmVal*ones(nu*4*M,1);
% y_ini = WarmVal*ones(nu*4*M,1);

% My own way
lambda_ini = max(lambda_ini,WarmVal);
y_ini = max(y_ini,WarmVal);

% % Some unclear method
% [m,~] = size(omega_r);
% y_ini = zeros(m,1);
% h_x = -OMEGA_L*delta_u_ini+omega_r;
% for i=1:m
% y_ini(i) = max([h_x(i),WarmVal]);
% end;


end