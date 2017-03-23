% This is a m-file to do Forward and back Substitution to solve Ly=b L'x=y;
% L comes form Cholesky Decomposition in order to solve Ax = b;
% The algorithm based on wikipedia - Triangular matrix
% DingYi. 2014.3.29


function x = lduEvaluate(L,D,U,b)
[m,n]=size(L);
y = zeros(m,1);
x = zeros(m,1);

L = L*D;
%Foward solve Ly=b
y(1) = b(1)/L(1,1);
for i=2:n
    y(i) = b(i)-L(i,1:i-1)*y(1:i-1);
    y(i) = y(i)/L(i,i);
end

%Backward solve Ux = y
x(n)=y(n)/U(n,n);
for j=1:n-1
    k = n-j;
    x(k) = y(k)-U(k,k+1:n)*x(k+1:n);
    x(k) = x(k)/U(k,k)';
end

end
        