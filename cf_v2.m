% This a function to do "Cholesky Factorization" A=L*L';
%   [L] = cf(A)
% Algorithm based on "Numerical Optimization", Algorithm 3.4
% DingYi
% 2014.03.29

function [L,d] = cf_v2(A)

[m,n] = size(A);
L = eye(m,n);
C = zeros(m,n);
d = zeros(m,1);
d(1,1) = A(1,1);
L(2:m,1) = A(2:m,1)/d(1,1);

for j=2:n
    C(j,j) = A(j,j)-d(1:j-1,1)'*(L(j,1:j-1).^2)';
    d(j,1) = C(j,j);
    for i=j+1:n
        C(i,j) = A(i,j)-d(1:j-1)'*(L(i,1:j-1).*L(j,1:j-1))';
        L(i,j) = C(i,j)/d(j,1);
    end
end
    
d = diag(d);

end