%测试自编的MPC函数是否可用
%2013.1.17
%注意，因为是初步的测试程序，因此在语句中有些是可以自动获取模型阶数的，有些是按照默认全量对象为2阶，增量对象为3阶写的
%因此要用到其他对象上时，需要再修改

%滚动时域控制的实现
clear;clc;
%建立对象模型（全量）
Ac = [0.5 1;0 0.5];
Bc = [0.5;1];
Cc = [1 0];
Dc = 0;

%预测控制参数配置
Np = 10;            %预测时域的选择，作阶跃测试，反向特性，纯滞后
Nc = 2;
R = 1;
Q = 1;

%augment是自定义的全量方程转增量方程的函数
[A_e,B_e,C_e] = augment(Ac,Bc,Cc,Dc);

%fphi是自定义的计算参数F和φ的函数，AA，BB为计算△U公式中的参数
[BarRs,Phi_Phi,Phi_F,Phi_R,F,Phi] = fphi(A_e,B_e,C_e,Np,Nc);
AA = (Phi_Phi+R*eye(Nc,Nc))\ Phi_R;
BB = (Phi_Phi+R*eye(Nc,Nc))\ Phi_F;

%滚动优化初始化
[n,n_in] = size(B_e);%确定状态维数
xm = [0;0];          %全量状态初始化
x_k = zeros(n,1);    %增量状态初始化
N_sim  = 50;         %仿真时域设定
r = ones(N_sim,1);   %目标曲线初始化
u_k = 0;            
y_k = 0;
u_k_1 = 0;
u_k_2 = 0;
delta_u1 = zeros(N_sim,1);
rr = ones(Np,N_sim);
delta_u2 = [0;0];

G = Phi'*(Q*eye(Np,Np))*Phi + R*eye(Nc,Nc);


%数组和矩阵的预先分配可以加快MATLAB运行速度
y1 = zeros(N_sim+1,1);
x1 = zeros(n,N_sim);
u1 = zeros(N_sim+1,1);

%tic;
%滚动优化
for kk = 1:N_sim;
    x_k = fankui(kk,x_k,y_k,u_k_1,u_k_2,A_e,B_e,C_e);
    c = (F*x_k-rr(:,kk))'*(Q*eye(Np,Np))*Phi;
    c = c';
    delta_u = quadprog(G,c,[1,0;-1,0;0,1;0,-1],0.8*[1;1;1;1],[],[],[],[],[0;0],'TolFun',1e-8);
    %[delta_u,~,~] = priduip(G,c,-[1,0;-1,0;0,1;0,-1],-0.8*[1;1;1;1],[0;0],[2;2;2;2],[0.5;0.1;0.1;0.1]);
    %[delta_u,~,~] = barrier(G,c,[1,0;0,-1],[0.3;0.3],[0;0]);
    delta_u = delta_u(1,1);
    %delta_u = youhua(AA,BB,r(kk),x_k);
    u_k = u_k + delta_u;
    delta_u1(kk) = delta_u;
    u1(kk+1,1) = u_k;
    xm = Ac*xm + Bc*u_k;  %全量状态更新
    y = Cc*xm;            %输出更新
    y1(kk+1,1) = y;
    y_k = y;              %将本周期内预测得到的y作为下一周期在线监测的y值
    x1(1,kk) = x_k(1,1);x1(2,kk) = x_k(2,1);x1(3,kk) = x_k(3,1);
    u_k_2 = u_k_1;
    u_k_1 = u_k;
end
%toc;

%画图
k=0:(N_sim-1);
kk=0:(N_sim);
%figure
%subplot(4,4,12);    plot(kk,y1,'LineWidth',2,'Color','r');title('y(k),AS,△U = 0.8');axis([0 15 0 1.5]);
%subplot(412);   plot(k,y2);  legend('y2');
%subplot(4,4,16);    stairs(kk,u1,'LineWidth',2);title('u(k),AS,△U = 0.8');axis([0 15 0 0.5]);
%subplot(414);   plot(k,u2);  legend('u2');
%xlabel('Sampling Instant');