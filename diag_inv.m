% This is m file to do diag inverse, because the inv() function in Matlab
% meets the problem of singular.
% Yi Ding 2014.4.26

function invMat = diag_inv(diag)
[m,n] = size(diag);
invMat = zeros(m,n);
for i=1:n
    invMat(i,i) = 1/diag(i,i);
end