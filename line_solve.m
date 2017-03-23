% This is a routine to solve linear equation using cholesky decomposition.
% Simply calling two other functions: cf() and luEvaluate()
% Act as the Matlab function "linsolve()" for Ax = b where A = A', A>0;
% DingYi. 2014.3.29


function x = line_solve(A,b)

L = cf(A);
%[L,D] = cf_v2(A);
x = luEvaluate(L,L',b);
%x = lduEvaluate(L,D,L',b);

end