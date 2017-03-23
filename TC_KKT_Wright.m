% One of the Termination Criteria: KKT residual (Wright)
% Based on the book Primal-Dual Interior-Point Methods, P226
% 2014.7.13
% Yi.

function flag = TC_KKT_Wright(rd,rp,b,c,mu,epsilon)

flag = 0;
temp1 = max(abs(rp))/(1+max(abs(b)));
temp2 = max(abs(rd))/(1+max(abs(c)));
 
if((temp1<=epsilon)&&(temp2<=epsilon)&&(mu<=epsilon))
    flag = 1;
end

end