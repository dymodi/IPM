% This a function to do "Cholesky Factorization" A=L*L';
%   [L] = cf(A)
% Algorithm based on "High Performance Cholesky and Symmetric Indefinite Factorizations with Applications" Algorithm 1
% DingYi
% 2013.12.05

function L = cf(A)

[m,n] = size(A);
L = zeros(m,n);

for j=1:n
    L(j,j) = sqrt(A(j,j));
    L(j+1:n,j)=A(j+1:n,j)/L(j,j);
    A(j+1:n,j+1:n)=A(j+1:n,j+1:n)-L(j+1:n,j)*L(j+1:n,j)';
end
