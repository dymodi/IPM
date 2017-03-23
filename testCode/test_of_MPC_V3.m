%改写原有MPC程序，优化结构，便于理解和代码移植
%2014.3.8
%本程序为测试程序，用来进行参数的设定和初始化，调用其他来完成MPC的算法
%暂时还是SISO，测试成功后。再找一个2x2的对象再来改写和测试。
%2014.4.15 用一个2x2的MIMO模型来测试 Test Successful!

clc;clear;
%变量维数确定
%暂时空缺，因为被控系统为SISO

%声明一些全局变量方便调用
global Ac Bc Cc A_e B_e C_e nu ny n_in;
global P M L;
global Phi F Q xm G GL;
global u_k_1 u_k_2;

global delta_U_p delta_U_n U_p U_n Y_p Y_n OMEGA_L;

global warm_y warm_lambda;

N_sim  = 50;                    %仿真时域设定

%建立对象模型（全量）
nu = 2;ny = 2;nx = 3;          %输入输出变量个数,其中nx为全量模型
Ac = [1,1,1;0,0,1;0,1,0];
Bc = [0,0;0,-1;1,0];
Cc = [0,0,1;0,0,1];
Dc = 0;

%预测控制参数配置
P = 5;
M = 2;
R = 1;
Q = 1;
L1 = zeros(nx,ny);
L2 = eye(ny,ny);
L=[L1;L2];  

%augment是自定义的全量方程转增量方程的函数，具体原理见文档
[A_e,B_e,C_e] = augment(Ac,Bc,Cc,Dc);

%fphi是自定义的计算参数F和φ的函数
[F,Phi] = fphi_v2(A_e,B_e,C_e,P,M);

%约束参数配置
II = eye(nu*M,nu*M);
B = eye(nu*M,nu*M);
for i = 1:nu*M
    for j = 1:nu*M
        for k = 1:M
            if(i==(j+(k-1)*nu))
                B(i,j)=1;
            end
        end
    end
end
%带被控变量的约束
%OMEGA_L = [II;-II;B;-B;Phi;-Phi];   
%不带被控变量的约束
OMEGA_L = [II;-II;B;-B];   
delta_U_p_ = [1;1];delta_U_n_ = [-1;-1];     %MIMO时这两个变量为向量
U_p_ = [0.4;0.4]; U_n_ = [-2;-2];            %MIMO时这两个变量为向量
Y_p_ = [2;2]; Y_n_ = [-2;-2];                %MIMO时这两个变量为向量
%将约束扩展到M（P）个时域
delta_U_p=[];delta_U_n=[];U_p=[];U_n=[];Y_p=[];Y_n=[];
for k = 1:M
    delta_U_p = [delta_U_p;delta_U_p_];
    delta_U_n = [delta_U_n;delta_U_n_];
    U_p = [U_p;U_p_];
    U_n = [U_n;U_n_];
end
for k = 1:P
    Y_p = [Y_p;Y_p_];
    Y_n = [Y_n;Y_n_];
end

%滚动优化初始化
[n,n_in] = size(B_e);           %确定状态维数
xm = zeros(nx,1);                     %全量状态初始化
x_k = zeros(n,1);               %增量状态初始化
r = ones(N_sim*ny,1);           %目标曲线初始化
u_k = zeros(nu,1);                        %控制变量初始化
y_k = zeros(ny,1);                        %被控变量初始化
delta_u = 0.1*ones(nu,1);             %用以给优化过程的x赋第一周期的初值；
u_k_1 = zeros(nu,1);                  %k-1时刻控制变量初始化
u_k_2 = zeros(nu,1);                  %k-2时刻控制变量初始化
delta_u2 = zeros(nu,1);
delta_u_M_in = zeros(nu*M,1);
delta_u_ini = 0.5*ones(nu*M,1);
y_ini = 0.5*ones(4*nu*M,1);
lambda_ini = 0.5*ones(4*nu*M,1);

%确定参考轨迹
rr = ones(ny*P,N_sim);

%G为二次规划求解中的参数：min 0.5*x'*G*x + c'*x   subject to:  A*x <= b
G = Phi'*(Q*eye(ny*P,ny*P))*Phi + R*eye(nu*M,nu*M);
GL = cf(G);

%记录各个变量在每个周期的值进行画图分析
delta_u_draw = zeros(N_sim,nu);
delta_u_uc_draw = zeros(N_sim,nu);
y_draw = zeros(N_sim+1,ny);
x_draw = zeros(n,N_sim);
u_draw = zeros(N_sim+1,nu);
Iter_rec = [];

%滚动优化
tic
for kk = 1:N_sim;
    r_k = rr(:,kk);        %对k时刻的参考轨迹进行更新
    [delta_u,delta_u_M_in,u_k,y_k,x_k,delta_u_ini,y_ini,lambda_ini] = mpc_v3(delta_u,delta_u_ini,u_k,y_k,x_k,r_k,delta_u_ini,y_ini,lambda_ini);%调用MPC在线算法进行计算
    %储存数据用于画图
    delta_u_draw(kk,:) = delta_u';
    u_draw(kk+1,:) = u_k';
    y_draw(kk+1,:) = y_k';          
    x_draw(1,kk) = x_k(1,1);x_draw(2,kk) = x_k(2,1);x_draw(3,kk) = x_k(3,1);x_draw(4,kk) = x_k(4,1);x_draw(5,kk) = x_k(5,1);  %这里默认的是有3个增广后的状态
    %Iter_rec(kk) = Iter;     %记录每周期优化算法的迭代次数
end
toc
%画图
figure;
subplot(2,2,1); plot(y_draw(:,1),'LineWidth',2,'Color','r'); title('y1(k)');%axis([0 15 0 1.5]);
subplot(2,2,2); plot(y_draw(:,2),'LineWidth',2,'Color','r'); title('y2(k)');%axis([0 15 0 1.5]);
subplot(2,2,3); stairs(u_draw(:,1),'LineWidth',2);           title('u1(k)');%axis([0 15 0 0.5]);
subplot(2,2,4); stairs(u_draw(:,2),'LineWidth',2);           title('u2(k)');%axis([0 15 0 0.5]);