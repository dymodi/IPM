% Algorithm "Modified Cholesky Factorization"
% L = mcf(A)
% According to the Numerical Optimization, 3.5 Modified Cholesky
% Factorization
% 
% Yi Ding 2014.5.6

function [L,D] = mcf(A)

[m,n] = size(A);
L = eye(m,n);
C = zeros(m,n);
d = zeros(m,1);
beta = 10;
sigma = 0.5;

C(1,1) = A(1,1);
d(1,1) = max([abs(C(1,1)),(C(1,1)/beta)^2,sigma]);
L(2:m,1) = A(2:m,1)/d(1,1);

for j = 2:n
    C(j,j) = A(j,j)-d(1:j-1)'*(L(j,1:j-1).^2)';
    d(j,1) = max([abs(C(j,j)),(max(abs(C(j:n,j)))/beta)^2,sigma]);
    for i=j+1:n
        C(i,j) = A(i,j)-d(1:j-1)'*(L(i,1:j-1).*L(j,1:j-1))';
        L(i,j) = C(i,j)/d(j,1);
    end
end

D = diag(d);

end