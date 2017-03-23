% This is a file to test the mimo system, we use a 2 input, 2output model.
% There are 3 state variables.

clc;clear;
A = [0,1,1;0,0,1;0,1,0];
B = [0,0;0,1;1,0];
C = [0,0,1;0,1,0];

u = ones(2,1);
y = zeros(2,1);
x = zeros(3,1);

Iter = 100;
x_k = zeros(3,Iter);
u_k = zeros(2,Iter);
y_k = zeros(2,Iter);

for k=1:Iter
    x_k(:,k) = x;
    u_k(:,k) = u;
    y_k(:,k) = y;
    x = A*x+B*u;
    y = C*x;
end

figure;
subplot(2,2,1);plot(u_k(1,:));
subplot(2,2,2);plot(u_k(2,:));
subplot(2,2,3);plot(y_k(1,:));
subplot(2,2,4);plot(y_k(2,:));
