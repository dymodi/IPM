%这个程序用来在Matlab上运行Model2的非线性的模型，以实现PIL
%相关的模型来自jMPC Toolbox中的nl_pend
%这个函数在Matlab_RS_V2的每次写周期中运行

function [y,xm2] = MPConDSP_mode2(u,xm2)

global Cc
global u_model2
x = xm2;

%Inputs
% u1 = Force applied to the cart [N]

%States
% x1 = Position of the cart [m]
% x2 = Velocity of the cart [m/s]
% x3 = Angle of the pendulum from vertical [rad]
% x4 = Angular velocity of the pendulum [rad/s]

%Assign Parameters
%[M,m,l,g] = param{:};


u_model2 = u;

%u = 0.1;

[T,x] = ode15s('model2ode',[0 0.05],x);

[row,col] = size(x);

xm2 = x(row,:)';

y = Cc*xm2;         


