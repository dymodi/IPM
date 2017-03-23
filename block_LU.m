% This is a function to do Block LU 
% Based on the Paper 《LU分解分块算法的研究与实现》
% Yi Ding
% 2014.4.25

function [L,U] = block_LU(A,mm,nn)
[m,n] = size(A);

A11 = A(1:mm,1:nn);
A12 = A(1:mm,nn+1:n);
A21 = A(mm+1:m,1:nn);
A22 = A(mm+1:m,nn+1:n);

[L11] = cf(A11);
U12 = A12/L11;
L21 = A21/L11';
[L22,U22,P] = lu(A22-L21*U12);

L = [L11,zeros(mm,n-nn);L21,P*L22];
U = [L11',U12;zeros(m-mm),U22];
end

