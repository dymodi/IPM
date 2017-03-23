%改写原有MPC程序，优化结构，便于理解和代码移植
%2014.3.8
%本程序为测试程序，用来进行参数的设定和初始化，调用其他来完成MPC的算法
%暂时还是SISO，测试成功后。再找一个2x2的对象再来改写和测试。
%2014.4.15 用一个2x2的MIMO模型来测试 Test Successful!
%2014.4.26 用文献《Auto-Code Generation for Fast Embedded Model Predictive
%Controllers》中的Model1来测试，比较效果
%2014.8.6 用文献《Auto-Code Generation for Fast Embedded Model Predictive
%Controllers》中的Model2来测试，比较效果
%2014.8.10用罗犇的Model3来测试

clc;clear;
%变量维数确定
%暂时空缺，因为被控系统为SISO

%声明一些全局变量方便调用
global Ac Bc Cc Cc1 A_e B_e C_e nu ny n_in ndec mc;
global P M L;
global Phi Phi1 Phi2 F F1 F2 Q QQ QQ2 xm G GL;
global u_k_1 u_k_2;
global QPTime TimeIter;
global delta_U_p delta_U_n U_p U_n Y_p Y_n OMEGA_L TotalIter;

global global_res;


N_sim  = 40;                    %仿真时域设定

%时间记录
QPTime = zeros(N_sim,1);
TotalIter = zeros(N_sim,1);
TimeIter = 1;

%建立对象模型（全量）
nu = 1;ny = 3;nx = 4;          %输入输出变量个数,其中nx为全量模型

Ac = [ 0.24    0     0.1787  0;
      -0.3722  1     0.2703  0;
      -0.9901  0     0.1389  0;
      -48.9354 64.1  2.3992  1];
%Bc = [-1.5429;-1.3648;-4.1106;11.4668];
Bc = [-1.2346;-1.4383;-4.4828;-1.7999];
Cc = [ 0     1      0  0;
       0     0      0  1;
      -128.2 128.2  0  0];
%Cc1 =[0 1 0 0] ;                %考虑约束时，把y1单独拉出来考虑，形成自己的Phi
Dc = 0;

%预测控制参数配置
P = 10;
M = 3;
R = 1;
Q = [1;1;1];          %还得修改，当MIMO时，R和Q都应该是向量
Q2 = 1;
%这里的L专门为积分对象进行修改
%L1 = eye(nx,ny);
L1 = zeros(nx,ny);
L2 = eye(ny,ny);
L=[L1;L2];  

%决策变量数和约束个数
ndec = nu*M;
%mc = 4*nu*M;        %不考虑输出变量的约束
mc = 4*nu*M+2*1*P; %考虑输出变量的约束

%augment是自定义的全量方程转增量方程的函数，具体原理见文档
[A_e,B_e,C_e] = augment(Ac,Bc,Cc,Dc);
%[A_e1,B_e1,C_e1] = augment(Ac,Bc,Cc1,Dc);   %用于单独对y1进行约束
C_e1 = C_e(1,:);
%C_e2 = C_e(2,:);

[F,Phi] = fphi_v2(A_e,B_e,C_e,P,M);
[F1,Phi1] = fphi_v2(A_e,B_e,C_e1,P,M);
%[F2,Phi2] = fphi_v2(A_e,B_e,C_e2,P,M);
%F1和Phi1的计算没问题，也理解了


%fphi是自定义的计算参数F和φ的函数
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
OMEGA_L = [II;-II;B;-B;Phi1;-Phi1];   
%不带被控变量的约束
%OMEGA_L = [II;-II;B;-B];   
delta_U_p_ = 0.524;delta_U_n_ = -0.524;     %MIMO时这两个变量为向量
U_p_ = 0.262; U_n_ = -0.262;                %MIMO时这两个变量为向量
Y_p_ = 0.349; Y_n_ = -0.349;       %MIMO时这两个变量为向量
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
delta_u_ini = 0.1*ones(nu*M,1);
y_ini = 0.5*ones(mc,1);
lambda_ini = 0.5*ones(mc,1);

%确定参考轨迹
rr1 = zeros(1,N_sim);
rr2 = 400*ones(1,N_sim);
rr3 = zeros(1,N_sim);
rr = []; 
for i=1:P
    rr = [rr;rr1;rr2;rr3];
end

%G为二次规划求解中的参数：min 0.5*x'*G*x + c'*x   subject to:  A*x <= b
QQ = [];
for k=1:P
    QQ = [QQ;Q];
    k=k+1;
end
QQ = diag(QQ);
QQ2 = Q2*ones(P,1);
QQ2 = diag(QQ2);
G = Phi'*QQ*Phi + R*eye(nu*M,nu*M);
GL = cf(G);

%记录各个变量在每个周期的值进行画图分析
delta_u_draw = zeros(N_sim,nu);
delta_u_uc_draw = zeros(N_sim,nu);
y_draw = zeros(N_sim+1,ny);
x_draw = zeros(n,N_sim);
u_draw = zeros(N_sim+1,nu);
Iter_rec = [];

%滚动优化
for kk = 1:N_sim;
    r_k = rr(:,kk);        %对k时刻的参考轨迹进行更新
    %tic
    [delta_u,delta_u_M_in,u_k,y_k,x_k,delta_u_ini,y_ini,lambda_ini] = mpc_v3_model3_full(delta_u,delta_u_ini,u_k,y_k,x_k,r_k,delta_u_ini,y_ini,lambda_ini);%调用MPC在线算法进行计算
    %QPTime(TimeIter,1) = toc;
    %TimeIter = TimeIter + 1;
    %储存数据用于画图
    delta_u_draw(kk,:) = delta_u';
    %delta_u_uc_draw(kk,:) = delta_u_uc';
    u_draw(kk+1,:) = u_k';
    y_draw(kk+1,:) = y_k';          
    %x_draw(1,kk) = x_k(1,1);x_draw(2,kk) = x_k(2,1);x_draw(3,kk) = x_k(3,1);%x_draw(4,kk) = x_k(4,1);x_draw(5,kk) = x_k(5,1);  %这里默认的是有3个增广后的状态
    %Iter_rec(kk) = Iter;     %记录每周期优化算法的迭代次数
end
%画图
u_draw = u_draw/3.1416*180;
figure;
subplot(4,1,1);plot(rr(1,:),'-.');hold on;plot(y_draw(:,1),'LineWidth',2,'Color','r');  title('Plant Output: yp1(k)'); %axis([0 N_sim -2 2]);
subplot(4,1,2);plot(rr(2,:),'-.');hold on;plot(y_draw(:,2),'LineWidth',2,'Color','r');  title('Plant Output: yp2(k)'); %axis([0 N_sim -2 2]);
subplot(4,1,3);plot(rr(3,:),'-.');hold on;plot(y_draw(:,3),'LineWidth',2,'Color','r');  title('Plant Output: yp3(k)'); %axis([0 N_sim -2 2]);
subplot(4,1,4);stairs(u_draw(:,1),'LineWidth',2);                                      title('Control Input:u(k)');  %axis([0 N_sim -5 5]);
hold on; plot(15*ones(N_sim,1),'--');hold on; plot(-15*ones(N_sim,1),'--');

figure;
subplot(2,1,1); stairs(1000*QPTime,'LineWidth',2);                                     title('QP Time /ms');%axis([0 15 0 0.5]);
subplot(2,1,2); stairs(TotalIter,'LineWidth',2);                                       title('QP Iterations');%axis([0 15 0 0.5]);

% figure;
% subplot(2,1,1);plot(delta_u_draw,'LineWidth',2);title('QP solved solution');
% subplot(2,1,2);plot(delta_u_uc_draw,'LineWidth',2);title('Uncontrained solution');
