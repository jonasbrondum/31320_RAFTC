%% 31320 Mandatory Assignment B
clear;
clc;
%% initial parameters
clear;
close all;
clc;
load('ECP_values.mat');
% Physical system parameters
J_1 = ECP_values(1);            % Disk 1 inertia kgm^2
J_2 = ECP_values(2);            % Disk 2 inertia kgm^2
J_3 = ECP_values(3);            % Disk 3 inertia kgm^2
k_1 = ECP_values(4);            % Shaft 1-2 stiffness Nm/rad
k_2 = ECP_values(5);            % Shaft 2-3 stiffness Nm/rad
b_1 = mean(ECP_values([6 7]));  % Disk 1 damping and friction Nms/rad
b_2 = mean(ECP_values([8 9]));  % Disk 2 damping and friction Nms/rad
b_3 = mean(ECP_values([10 11]));% Disk 3 damping and friction Nms/rad
T_Cp = ECP_values(12);          % Disk 1 Coulomb friction in positive direction
T_Cm = ECP_values(13);          % Disk 1 Coulomb friction in negative direction
atan_scale = 100;               % Sign approximation factor
w_th = 0.75;                    % Threshold angular velocity rad/s

%% State space representation
% x = [theta1, omega1, theta2, omega2, theta3, omega3] = [x1 .. x6]
syms x1 x2 x3 x4 x5 x6 d1 u1 u2
f1 = x2;
f2 = 1/J_1*((u1) - b_1*x2 - k_1*(x1 - x3) - d1);
f3 = x4;
f4 = 1/J_2*((u2) - b_2*x4 - k_1*(x3 - x1) - k_1*(x3 - x5));
f5 = x6;
f6 = 1/J_3*(-b_3*x6 - k_2*(x5 - x3));
f = [f1; f2; f3; f4; f5; f6];
g1 = x1;
g2 = x3;
g3 = x5;
g = [g1; g2; g3];
x = [x1; x2; x3; x4; x5; x6];
u = [u1; u2];
d = [d1];
A = double(jacobian(f, x));
B = double(jacobian(f, u));
C = double(jacobian(g, x));
D = double(jacobian(g, u));
E_x = double(jacobian(f, d));
E_y = double(jacobian(g, d));


%% SISO system from input to top disk
[num, den] = ss2tf(A,B,C,D,1);

% Top disk is y3 which is theta_3:
% SISO system from state-space
G = tf(num(3,:),den);

%% Sampling period

T_s = 0.004;                    % Sampling period
%% Redefinition of weights for Simulink simulation

% First we need to make some weights relevant as to the design
% requirements:

% We want gamma ~= 1
% And we don't want too high gains
% We want (almost) integral action meaning (almost) no steady state error
% We don't want the actuator to go into saturation for a step-response

% Example from book
% Uses the Robust Control toolbox
M=2; wb=5; A=1.e-3; % Hvad betyder disse og hvor i bogen kommer de fra?
%Ligning (2.113)
%M=10; wb=3; A=1.e-4; % Hvad betyder disse og hvor i bogen kommer de fra?
W1 = tf([1/M wb], [1 wb*A]);


[A,B,C,D] = tf2ss(W1.Numerator{1},W1.Denominator{1});

W1 = ss(A,B,C,D);
%W1= makeweight(20,35,0.1)

% W2=tf(0.008,1)
% 
% [A,B,C,D] = tf2ss(W2.Numerator{1},W2.Denominator{1});
% 
% W2 = ss(A,B,C,D);

W2 = makeweight(0.001,40,20) %setup 1

[K,CL,gamma] = mixsyn(G,W1,W2,[],1); %Last argument makes the function try to force gamma to 1
% [K,CL,gamma] = mixsyn(G,W1,[],[]);
gamma

% First, compare the resulting sensitivity S and complementary sensitivity 
% T to the corresponding weighting functions W1 and W3. 

L = G*K;
I = eye(size(L));
SenFun = feedback(I,L); 
T= I-SenFun;

figure;
sigma(SenFun,'b',W1,'b--',T,'r',W2,'r--',{0.01,1000})
legend('S','W1','T','W2','FontName','times','Interpreter','latex')
title('Sensitivity and weight functions','Interpreter','latex')
grid on
saveas(gcf,'figures/Q9_sensitivity_weight_functions.svg')

figure;
sigma(L,'b',W1,'r--',1/W2,'g--',{0.01,1000})
legend('L','W1','1/W2','Interpreter','latex')
title('Open-loop and weight functions','Interpreter','latex')
saveas(gcf,'figures/Q9_open_loop_weight_functions.svg')
close all

%Check phasemargin
Marg = allmargin(G*K)

%% Test of implementation:
% Augmented plant for synthesis
%                                       +------+
%            +--------------------------|W_1(s)|----> z1
%            |                          +------+
%            |                        +------+
%            |                 +------|W_2(s)|------> z2
%            |                 |      +------+
%         +  | e +------+  u(s)| +------+             
% ref ---->O-----| K(s) |--------| G(s) |-----------> y
%         -|     +------+        +------+       |         
%          |                                    |
%          +------------------------------------+
%                  

G.InputName = 'ug';
G.OutputName= 'yg';

W1.InputName = 'uW1';
W1.OutputName = 'yW1';

W2.InputName = 'uW2';
W2.OutputName = 'yW2';

K.InputName = 'uK';
K.OutputName = 'yK';

blksys = append(G,W1,W2,K);

% specify connections
S1 = sumblk('ug = yK',1);
S2 = sumblk('y = yg',1);
S3 = sumblk('uW1 = e',1);
S4 = sumblk('uK = e',1);
S5 = sumblk('z1 = yW1',1); 
S6 = sumblk('e = r - yg',1);
S7 = sumblk('uW2 = yK', 1);
S8 = sumblk('z2 = yW2',1);
S9 = sumblk('u = yK',1);
S = append(S1,S2,S3,S4,S5,S6,S7,S8,S9);
 


 % connect to obtain generalised plant
% P    = [Wp , Wp*G, Wp*G;
%         0  , 0   , Wu;
%         I  , G   , G];
Pcl  = connect(blksys,S,{'r'},{'y','z1','z2','e','u'});

Gcl = feedback(K*G,1);
% % validate: 
% figure ,bode(Pcl('y','r'), (G*K)/(1+G*K),'--r',{0.1,1000})
% disp ('validate that G*K/(1+G*K) == transfer function from ref to y')
% figure ,bode(Pcl('z1','r'), W1/(1+G*K),'--r',{0.1,1000})
% disp ('validate that W1/(1+G*K) == transfer function from ref to z1')
% figure ,bode(Pcl('z2','r'), W2/(1+G*K),'--r',{0.1,1000})
% disp ('validate that W2/(1+G*K) == transfer function from ref to z2')
% figure ,bode(Pcl('e','r'), 1/(1+G*K),'--r',{0.1,1000})
% disp ('validate that 1/(1+G*K) == transfer function from ref to e')

% Controller assumes positive feedback

hinfnorm(Pcl('y','r'))
hinfnorm(Pcl('z1','r'))
hinfnorm(Pcl('z2','r'))
hinfnorm(Pcl('e','r'))
 
% figure;
% step((Pcl))
% figure;

% step(c2d(Pcl,T_s))


%% Simulation in Simulink
%Simulink outputs very different results from Henriks method. This is
%probably due to his system being continous, while the simulink model is
%discrete


[numK, denK]=ss2tf(K.A, K.B, K.C, K.D,1)

[numW1, denW1]=ss2tf(W1.A, W1.B, W1.C, W1.D,1)


[numW2, denW2]=ss2tf(W2.A, W2.B, W2.C, W2.D,1)



x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.0004;                    % Sampling period
