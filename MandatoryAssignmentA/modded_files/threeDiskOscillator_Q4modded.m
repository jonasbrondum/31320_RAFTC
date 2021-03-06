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

% The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;                    % Sampling period
sigma_meas = 0.0093*eye(3);     % Measurements covariance matrix


% simTime = 1;

%% SIMULINK residuals continuous time implementation:
s = tf('s')



%% State space representation
% x = [theta1, omega1, theta2, omega2, theta3, omega3] = [x1 .. x6]
syms x1 x2 x3 x4 x5 x6 d1 u1 u2 fa1 fa2 fa3 fa4 fa5
f1 = x2;
f2 = 1/J_1*((u1 + fa1) - b_1*x2 - k_1*(x1 - x3) - d1);
f3 = x4;
f4 = 1/J_2*((u2 + fa2) - b_2*x4 - k_1*(x3 - x1) - k_1*(x3 - x5));
f5 = x6;
f6 = 1/J_3*(-b_3*x6 - k_2*(x5 - x3));
f = [f1; f2; f3; f4; f5; f6];
g1 = x1 + fa3;
g2 = x3 + fa4;
g3 = x5 + fa5;
g = [g1; g2; g3];
x = [x1; x2; x3; x4; x5; x6];
u = [u1; u2];
d = [d1];
faults = [fa1; fa2; fa3; fa4; fa5];
A = double(jacobian(f, x));
B = double(jacobian(f, u));
C = double(jacobian(g, x));
D = double(jacobian(g, u));
E_x = double(jacobian(f, d));
E_y = double(jacobian(g, d));
F_x = double(jacobian(f, faults));
F_y = double(jacobian(g, faults));

% Discrete time
sys_d = c2d(ss(A, B, C, D), T_s);
F_d = sys_d.A;
G_d = sys_d.B;

% State-feedback LQR design
Q_c = diag([2 0 2 0 2.5 0.0024]);
R_c = diag([10 10]);
K_c = [];

% Scaling of reference
C_ref = [];

% Kalman filter with friction estimation - DO NOT MODIFY
F_aug = [F_d G_d(:,1);zeros(1,6) 1];
G_aug = [G_d;0 0];
C_aug = [C zeros(3,1)];
% Kalman gain
L_aug = dlqe(F_aug,eye(7),C_aug,1e-3*eye(7),deg2rad(0.0056)*eye(3));
L_o = L_aug(1:6,:);
L_d = L_aug(7,:);

%% H_rf(s)
% This transfer function from fault to residuals is first of all the
% product: V_ry(s)*H_yf(s) = H_rf(s)
% H_yf(s) = C*(s*I - A)^(-1)*F_x + F_y
% Also F(s) = [V_ry(s) V_ru(s)]
% F = null(H')'
% H = [H_yu(s) H_yd(s);
%      I       0    ]
% H_yu(s) = C*(s*I - A)^(-1)*B + D
% H_yd(s) = C*(s*I - A)^(-1)*E_x + E_y
syms s
I = eye(size(A));
H_yf = C*(s*I - A)^(-1)*F_x + F_y;
H_yu = C*(s*I - A)^(-1)*B + D;
H_yd = C*(s*I - A)^(-1)*E_x + E_y;
H = [H_yu H_yd;
     eye(2) zeros(2,1)];
F = (simplify(null(H')'));
V_ry = F(:,1:3)
H_rf = vpa(simplify(V_ry*H_yf),4)
%% Residual filter design

%% Strong and weak detectability
% H_rf = tf(0);
% s = tf('s');

% Weak detectabilities:
% The ith fault is weakly detectable if and only if:
% rank([ H_yd H^i_yf]) > rank(H_yd)
dimH_yf = size(H_yf)
weakly_detectable = zeros(dimH_yf(2));
for i = 1:dimH_yf(2)
    if rank([H_yd H_yf(:,i)]) > rank(H_yd)
        txt = sprintf('Fault %d is weakly detectable',i);
        disp(txt)
        weakly_detectable(i) = 1;
    end
end

% Now we check also if these are strongly detectable:
% Strong detectability:
% ith fault is strongly detectable if F(s)*[H^i_yf; 0]|s=0 != 0
% so if that product with s evaluated at 0 is not zero then the fault is
% also strongly detectable, aka it has a steady state gain different from
% zero.
added_zero_vector = zeros(2,1);
s = 0;
strongly_detectable = zeros(dimH_yf(2));
for i = 1:dimH_yf(2)
    if eval(F*[H_yf(:,i); added_zero_vector]) ~= added_zero_vector
        txt = sprintf('Fault %d is strongly detectable',i);
        disp(txt)
        strongly_detectable(i) = 1;
    end
end
%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
h = 0;                  % Put the threshold from GLR here

%% Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
va_eig_d = [];  % Discrete time eigenvalues
va_eig = log(va_eig_d)/T_s;     % Continuous time eigenvalues
% Then discretise your VA

B_change = [1 0;0 0];

%% Simulation for sensor fault (f_u = 0)
simTime = 45;                   % Simulation duration in seconds
% f_u_time = 25;                  % Actuator fault occurence time
% detect_time = f_u_time + 3.75;
% f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
% u_fault = 0;                    % Disable VA meachanism
% f_m_time = 8.5;                 % Sensor fault occurence time
% sim('threeDiskOscillatorRig');

%% Simulation for actuator fault (f_m = 0)
% f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
% u_fault = 1;                    % Enable VA meachanism
% f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
% sim('threeDiskOscillatorRig_solution');
%  
%% Plot settings
% set(0,'DefaultTextInterpreter','latex');
% set(0,'DefaultAxesFontSize',20);
% set(0,'DefaultLineLineWidth', 2);

