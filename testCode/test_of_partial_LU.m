%This is a test file to test whether partial factorlzation can be
%implemented
%2014.4.25
%Yi Ding

clc;clear;

G = [38.1,-26,26,-20;-26,20.1,-20,14;26,-20,20.1,-14;-20,14,-14,12.1];
A = [-1,1,0,0;0,-1,1,0;0,0,-1,1;0,0,0,-1];
lambda_y = [0.1,0,0,0;0,0.2,0,0;0,0,0.3,0;0,0,0,0.4];

F = [G,A';A,-lambda_y];
[L,U,P] = lu(F);
[Lb,Ub] = block_LU(F,4,4);