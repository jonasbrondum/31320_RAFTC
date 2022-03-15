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

%% Filtering

s = tf('s')

sys = -s^2-(b_2*s-k_2-k_1)/J_2
zeta=0.7;
% Originally we chose omega = 10/4
omega=10; % From experiment 28/02 -> omega = 20 was a good value


lowpass1=omega^2/(s^2+2*zeta*omega*s+omega^2)


%% Discretize

s = tf('s')


F=[ 0 1/J_2 k_1/J_2 (-s^2*J_2-b_2*s-k_1-k_2)/J_2 k_2/J_2;
    0 0 0 k_2/J_3 (-J_3*s^2-b_3*s-k_2)/J_3];
    
%Q = [lowpass1 0;
%    0 lowpass1]

RG1 = J_2*c2d(lowpass1*F(1,:),T_s,'tustin')
RG2 = J_3*c2d(lowpass1*F(2,:),T_s,'tustin')

[num1 den1] = tfdata(RG1);
[num2 den2] = tfdata(RG2);

%residual 3
RG3=c2d(lowpass1*F(1,4)/F(2,4),T_s,'tustin')
[num3 den3] = tfdata(RG3);

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
% All is denominated: [u1 u2 y1 y2 y3] -> [fa1 fa2 fa3 fa4 fa5]

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
% s = 0;
strongly_detectable = zeros(dimH_yf(2));
for i = 1:dimH_yf(2)
    disp(eval(subs(F*[H_yf(:,i); added_zero_vector],s,0)))
    if eval(subs(F*[H_yf(:,i); added_zero_vector],s,0)) ~= added_zero_vector
        txt = sprintf('Fault %d is strongly detectable',i);
        disp(txt)
        strongly_detectable(i) = 1;
    end
end

% D)
% If the disturbance was known then at least there would be another
% constraint available as a residual. This could result in better
% detectability of the faults both in terms of u1, but also in terms of
% weak versus strong detectability.
%% Variance of residuals:
% First we convert transfer function to state-space
sysss = ss(d2c(RG2))
sigma_y = 0.0093;
Ar = sysss.A;
Br = sysss.B;
Cr = sysss.C;
Dr = sysss.D;
% The covariance of the "inputs" [u1 u2 y1 y2 y3] are given:
Q_wr = sigma_y^2*eye(5); % Should be squared

syms q1 q2 q3 q4
% Now we do Lyapunov's equation to find Q_xr
% Q_xr = [q1 q2;
%         q3 q4];
% eq_lyap = 0*eye(length(Ar)) == Ar*Q_xr + Q_xr*(Ar') + (Br*Q_wr*(Br'));

% [q1, q2, q3, q4] = vpasolve(eq_lyap, [q1, q2, q3, q4])
% Q_test = Br*Q_wr*(Br');
% 
% Q_xr = dlyap(Ar, Q_test);
% 
% 
% % Q_xr = [q1 q2;
% %         q3 q4];
% 
% Q_yr = Cr*Q_xr*Cr' + Dr*Q_wr*Dr';
%     
% % Since this is sigma_yr^2 we have to take the square root:
% sigma_yr = double(sqrt(Q_yr))

Q_xr = lyap(Ar, Br*Q_wr*Br');
Q_yr = Cr*Q_xr*Cr' + Dr*Q_wr*Dr';

sigma_yr = sqrt(Q_yr)
%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])

% False alarm probability:
P_F = 0.0001;

% Probability of missed detection
P_M = 0.01;

% Probability of detection
P_D = 1 - P_M;


h = chi2inv(1 - P_F, 1)/2;                  % Put the threshold from GLR here

% We want to make sure that the threshold P_F = 1 - chi2cdf(2*h, 1) and
% -> h = chi2inv(1 - P_F, 1)/2
% For the window size it's slightly more of a hassle:
% P_D = 1 - P_M
syms X M;     % zz is the integrant variable (X), gg is lambda in symbolic

mu_1 = dcgain(RG2(4))*f_m(2);
mu_0 = 0;

lambda = (M*(mu_1 - mu_0)^2)/(sigma_yr^2);

p_x2 = 0.5*(X/lambda)^(-0.25)*exp(-(X + lambda)/2)*besseli(-0.5, sqrt(lambda*X));
integral = int(p_x2, X, 2*h, inf);

M = ceil(double(vpasolve([integral == P_D], [M], [0 200])))


r_window=zeros(M,1);
 

% Density function expression
% pd_zz = ( (1/2)*(zz/gg)^(-0.25) )* exp(-(zz + gg)/2)*besseli(-0.5, sqrt(gg*zz));
% p_zz = int(pd_zz, zz,2*h, Inf);  % FILL IN - Integrate over the probability space
% eq_1 = p_zz == P_D;  % FILL IN - Equation to be solved
% lambda = double(vpasolve(eq_1, gg)); %1.2625
% lambda=1.2625;
%% Test lambda
clc;
% lambda = 0:0.1:100;
% pTest = ncx2cdf(2*h, 1, lambda);
% pTest2 = 1 - ncx2cdf(2*h, 1, lambda);
for i = 1:100
    pTest = ncx2cdf(2*h, 1, i*(mu_1^2)/(sigma_yr^2));
    if pTest < P_M
        M2 = i
        break
    end
end
% hold on
% plot(lambda, pTest)
% plot(lambda, pTest2)
% disp(lambda)
% M = lambda*sigma_yr^2/((mu_1 - mu_0)^2)
% plot(lambda, M)

%% Q6

load('Experiment\Second_run.mat');

figure
hold on
plot(g_GLR,'LineWidth',2)
plot(h_GLR,'LineWidth',2)
legend('$\mathbf{g}(h)$', '$\mathbf{h}$','FontSize',16,'Interpreter','latex');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{g}(h)$','FontName','times','FontSize',16,'Interpreter','latex')
title('GLR test during sensor fault');
hold off

%% 

figure
hold on
plot(H_GLR,'LineWidth',2)
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{H}(h)$','FontName','times','FontSize',16,'Interpreter','latex')
title('Hypothethis during GLR test');
hold off



%% Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
va_eig_d = [];  % Discrete time eigenvalues
va_eig = log(va_eig_d)/T_s;     % Continuous time eigenvalues
% Then discretise your VA

B_change = [1 0;0 0];

%% DLQR


Q_c = [2, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 2, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 2.5, 0;
       0, 0, 0, 0, 0, 0.0024];
R_c =[10, 0;
      0, 10];

[K_c, S, CLP]=dlqr(F_d, G_d, Q_c, R_c);

C_3 = [0 0 0 0 1 0];

% State-feedback LQR design
% Q_c = diag([2 0 2 0 2.5 0.0024]);
% R_c = diag([10 10]);
% K_c = [];

% Scaling of reference
% C_ref = [];

% Kalman filter with friction estimation - DO NOT MODIFY
F_aug = [F_d G_d(:,1);zeros(1,6) 1];
G_aug = [G_d;0 0];
C_aug = [C zeros(3,1)];
% Kalman gain
L_aug = dlqe(F_aug,eye(7),C_aug,1e-3*eye(7),deg2rad(0.0056)*eye(3));
L_o = L_aug(1:6,:);
L_d = L_aug(7,:);
%% 


C_ref=pinv(C_3*( eye(6)-F_d + G_d*K_c )^(-1) * G_d*K_c);


%% Q8 - Virtual Actuator
% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
B_f = B(:,1)

% 1. Can there be perfect static matching?
if rank(B_f) == rank([B B_f])
    disp('Perfect static matching is possible')
else
    disp('Imperfect static matching is possible')
end

% There is imperfect static matching as the condition above isn't fulfilled
% This means we will have to sort to a dynamic reconfiguration.

% 2. Simulation of fault with DLQR - no fault accomodation
% Then design discret virtual actuator:
G_f = G_d(:, 1);
if rank(G_f) == rank([G_d G_f])
    disp('Perfect static matching is possible')
else
    disp('Imperfect static matching is possible')
end

% We check controllability
if max(size(F_d)) == rank(ctrb(F_d, G_f))
    disp('Faulty system is controllable');
else
    disp('Faulty system is not controllable');
end
% It is controllable

% Continuous time
va_eig = log(eig(F_d - G_d*K_c))/T_s;
M_c = place(A, B_f, va_eig);
A_D = A - B_f*M_c;
N_D = pinv(B_f)*B;
B_D = B - B_f*N_D;
C_D = C;

% Discrete time
va_eig_d = exp(va_eig*T_s);
M_d = place(F_d, G_f, va_eig_d);
F_D = F_d - G_f*M_d;
N_D_d = pinv(G_f)*G_d;
G_D = G_d - G_f*N_D_d;
C_D = C;


% 3. Implement and simulate with VA

% Block diagram is implemented in threeDiskOscillatorRig_Q6.slx



% sim('threeDiskOscillatorRig_solution');

%% Virtual sensor
dimA = max(size(A));
E = zeros(dimA, 1);
E(2) = 1;
A_aug = [A E;
         zeros(1, dimA) 0];
B_aug = [B; zeros(1, 2)];
C_aug = [C zeros(3, 1)];
rowC = max(size(x));
colC = min(size(C));
C_f = C;
C_f(2,:) = zeros(1, rowC);
C_f_aug = [C_f zeros(colC, 1)];

if rank(C_f_aug) == rank([C_aug; C_f_aug])
    disp('Perfect static matching for sensor fault');
else
    disp('Imperfect static matching for sensor fault');
end

if max(size(A_aug)) == rank(obsv(A_aug, C_f_aug)')
    disp('Faulty system is observable');
else
    disp('Faulty system is not observable');
end
% Continuous time
vs_eig = [-20 -41 -52 -45 -23 -25 -25];
L_V = place(A_aug', C_f_aug', vs_eig)';
A_V = A_aug - L_V*C_f_aug;
B_V = B_aug;
P_V = C_aug*pinv(C_f_aug);
C_V = C_aug - P_V*C_f_aug;

% Discrete time
F_aug = [F_d E;
         zeros(1, dimA) 1];
G_aug = [G_d; zeros(1, 2)];
C_f_aug = C_f_aug;
vs_eig_d = exp(vs_eig*T_s);
L_V_d = place(F_aug', C_f_aug', vs_eig_d)';
F_V =  F_aug - L_V_d*C_f_aug;
G_V = G_aug;
P_V_d = C_aug*pinv(C_f_aug);
C_V_d = C_aug - P_V_d*C_f_aug;

%save('Residual_generatorsQ2.mat');

%sim('threeDiskOscillatorRig');

%% Simulation for sensor fault (f_u = 0)


f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m_time = 8.5;                 % Sensor fault occurence time


%% Plot settings
% set(0,'DefaultTextInterpreter','latex');
% set(0,'DefaultAxesFontSize',20);
% set(0,'DefaultLineLineWidth', 2);

%% Simulating without virtual actuator
% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
detect_time=45;

figure
hold on
plot(theta_ref,'LineWidth',2)
plot(y_contr,'LineWidth',2)
legend('$\mathbf{\theta_{ref}}$', '$\mathbf{\theta_{measured}}$','FontSize',16,'Interpreter','latex');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{\theta}(t) [rad]$','FontName','times','FontSize',16,'Interpreter','latex')
title('LQR controller actuator fault');
hold off
%% With virtual actuator
% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
detect_time = f_u_time + 3.75;

figure
hold on
plot(theta_ref,'LineWidth',2)
plot(y_contr,'LineWidth',2)
legend('$\mathbf{\theta_{ref}}$', '$\mathbf{\theta_{measured}}$','FontSize',16,'Interpreter','latex');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{\theta}(t) [rad]$','FontName','times','FontSize',16,'Interpreter','latex')
title('Virtual actuator performance');
hold off


%% Residuals when actuator fault


figure
hold on
plot(RG1,'LineWidth',2)
plot(RG2,'LineWidth',2)
plot(RG3,'LineWidth',2)
legend('RG1', 'RG2', 'RG3');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{r}(t)$','FontName','times','FontSize',16,'Interpreter','latex')
hold off
