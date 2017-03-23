%改写原有MPC程序，优化结构，便于理解和代码移植
%2014.3.8
%本程序为测试程序，用来进行参数的设定和初始化，调用其他来完成MPC的算法
%暂时还是SISO，测试成功后。再找一个2x2的对象再来改写和测试。

clc;clear;
%变量维数确定
%暂时空缺，因为被控系统为SISO

%声明一些全局变量方便调用
global Ac Bc Cc A_e B_e C_e nu ny n_in;
global P M L;
global Phi F Q xm G;
global u_k_1 u_k_2;

global delta_U_p delta_U_n U_p U_n Y_p Y_n OMEGA_L;

%建立对象模型（全量）
nu = 1;ny = 1;nx = 2;          %输入输出变量个数,其中nx为全量模型
Ac = [0.5 1;0 0.5];
Bc = [0.5;1];
Cc = [1 0];
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
[BarRs,Phi_Phi,Phi_F,Phi_R,F,Phi] = fphi(A_e,B_e,C_e,P,M);

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
OMEGA_L = [II;-II;B;-B;Phi;-Phi];   
%不带被控变量的约束
%OMEGA_L = [II;-II;B;-B];   
delta_U_p_ = 0.6;delta_U_n_ = -0.2;    %MIMO时这两个变量为向量
U_p_ = 0.6; U_n_ = -1;                 %MIMO时这两个变量为向量
Y_p_ = 1.02; Y_n_ = -1;                %MIMO时这两个变量为向量
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
xm = [0;0];                     %全量状态初始化
x_k = zeros(n,1);               %增量状态初始化
N_sim  = 50;                    %仿真时域设定
r = ones(N_sim,1);              %目标曲线初始化
u_k = 0;                        %控制变量初始化
y_k = 0;                        %被控变量初始化
delta_u = 0.1;                  %用以给优化过程的x赋第一周期的初值；
u_k_1 = 0;                      %k-1时刻控制变量初始化
u_k_2 = 0;                      %k-2时刻控制变量初始化
delta_u2 = [0;0];

%确定参考轨迹
rr = ones(P,N_sim);

%G为二次规划求解中的参数：min 0.5*x'*G*x + c'*x   subject to:  A*x <= b
G = Phi'*(Q*eye(P,P))*Phi + R*eye(M,M);

%记录各个变量在每个周期的值进行画图分析
delta_u_draw = zeros(N_sim,1);
y_draw = zeros(N_sim+1,1);
x_draw = zeros(n,N_sim);
u_draw = zeros(N_sim+1,1);
Iter_rec = [];

%滚动优化
for kk = 1:N_sim;
    r_k = rr(:,kk);        %对k时刻的参考轨迹进行更新
    [delta_u,u_k,y_k,x_k,Iter] = mpc_v2(delta_u,u_k,y_k,x_k,r_k);%调用MPC在线算法进行计算
    %储存数据用于画图
    delta_u_draw(kk) = delta_u;
    u_draw(kk+1,1) = u_k;
    y_draw(kk+1,1) = y_k;          
    x_draw(1,kk) = x_k(1,1);x_draw(2,kk) = x_k(2,1);x_draw(3,kk) = x_k(3,1);  %这里默认的是有3个增广后的状态
    Iter_rec(kk) = Iter;     %记录每周期优化算法的迭代次数
end

%画图
figure;
subplot(2,1,1); plot(y_draw,'LineWidth',2,'Color','r'); title('y(k)');%axis([0 15 0 1.5]);
subplot(2,1,2); stairs(u_draw,'LineWidth',2);           title('u(k)');%axis([0 15 0 0.5]);
