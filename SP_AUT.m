% This is one of the ways to choose a starting point -- AUT's method.
% Based on the AUT's code: quad_wright
% 2014.7.13
% Yi

function [delta_u_ini,y_ini,lambda_ini] = SP_AUT(G,c,A,b,WARMVAL)

ndec = length(c); 
mc = length(b);

pmax = max(max(abs([G c; A b])));
if(pmax > 1)
    WARMVAL = sqrt(pmax);
end

delta_u_ini = zeros(ndec,1);
y_ini = WARMVAL*ones(mc,1);
lambda_ini = WARMVAL*ones(mc,1);

end