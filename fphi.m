%fphi是自定义的计算参数F和φ的函数
function [F,Phi] = fphi(A_e,B_e,C_e,Np,Nc)
[m1,n1] = size(A_e);
[m2,n2] = size(B_e);
[m3,n3] = size(C_e);

%计算F阵与φ阵
%F = zeros(Np*m3,n1);
F = [];
for i=1:Np
    F = [F;C_e * A_e^i];
end
v = zeros(Np,1);
for(i=1:Np)
    v(i,1) = C_e*A_e^(i-1)*B_e;
end
Phi = zeros(Np,Nc);
Phi(:,1) = v;
for(i=2:Nc)
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)];
end
% BarRs=ones(Np,1);
% Phi_Phi= Phi'*Phi;
% Phi_F= Phi'*F;
% Phi_R=Phi'*BarRs;