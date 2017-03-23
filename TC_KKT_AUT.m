% One of the Termination Criteria: KKT residual (AUT)
% Based on the quad_wright code.
% 2014.7.14
% Yi.

function flag = TC_KKT_AUT(rd,rp,p_max,mu,epsilon)

flag = 0;

temp1 = max([abs(rp);abs(rd)]);
phi = temp1/p_max;
 
if((phi<=epsilon)&&(mu<=epsilon))
    flag = 1;
end

end