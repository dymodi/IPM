%优化环节，即计算增量输出
%delta_u的计算基于【书籍】《Model Predictive Control System Design and Implementation
%Using MATLAB》P13, 1.25
function [delta_u] = youhua(AA,BB,r_k,x_k)
delta_U = AA * r_k - BB * x_k;
delta_u = delta_U(1,1);
