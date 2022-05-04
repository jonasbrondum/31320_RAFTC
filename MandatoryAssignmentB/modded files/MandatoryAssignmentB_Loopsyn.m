%% 31320 - Mandatory Assignment B
% Question 13 - with multiplicative uncertainty

load('ECP_values.mat');
% Physical system parameters
J_1 = 0.0325;  %CHANGED
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
G_p = tf(num(3,:),den);
%% Redefinition of weights for Simulink simulation

load('G.mat');
s = tf ('s')

Gd = tf(10,[1 0]);

G_p=6.498e07/(s^6 + 0.276*s^5 + 4520*s^5 + 830.7 * s^4 + 2.31e06*s^3 +2.1e05 * s^2 +21e05*s); %removed last element which was miniscule

alpha = 0;
[K_mix,CL_mix,gamma_mix,info_mix] = loopsyn(G_p,Gd);
gamma_mix


