% One of the Termination Criteria: Himmelblau
% Based on the Jiaona Wan's Doctor thesis
% 2014.7.14
% Yi.

function [flag,fx_old,x_old] = TC_Himmelblau(G,c,A,b,fx_old,x_old,x,epsilon)
%global global_res;

flag = 0;

% Himmelblau Termination Criteria
fx = 0.5*x'*G*x + c'*x;
res1 = abs(fx_old - fx);
res2 = abs(x_old - x);
cons_obey = min(A*x-b);
res = max([res1;res2]);
%global_res(Iter+1) = res;
fx_old = fx;
x_old = x;
% Check residuel
if((cons_obey>=-epsilon)&&(res<epsilon))
    flag = 1;
end

end