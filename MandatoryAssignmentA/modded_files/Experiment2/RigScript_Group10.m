clear all;
load('rig_init.mat');

%% Filtering

s = tf('s')

sys = -s^2-(b_2*s-k_2-k_1)/J_2
zeta=0.7;
omega=sqrt(k_2/J_3)
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


num12=cell2mat(num1(2));
num13=cell2mat(num1(3));
num14=cell2mat(num1(4));
num15=cell2mat(num1(5));



den12=cell2mat(den1(2));
den13=cell2mat(den1(3));
den14=cell2mat(den1(4));
den15=cell2mat(den1(5));



num24=cell2mat(num2(4));
num25=cell2mat(num2(5));
den24=cell2mat(den2(4));
den25=cell2mat(den2(5));



%% Variance of residuals:
% First we convert transfer function to state-space
sysss = ss(RG2)
sigma_y = 0.0093;
Ar = sysss.A;
Br = sysss.B;
Cr = sysss.C;
Dr = sysss.D;
% The covariance of the "inputs" [u1 u2 y1 y2 y3] are given:
Q_wr = sigma_y^2*eye(5); % Should this be squared or not???? <-----
% Q_wr(1:2,1:2) = 0;

% syms q1 q2 q3 q4
% % Now we do Lyapunov's equation to find Q_xr
% Q_xr = [q1 q2;
%         q3 q4];
% eq_lyap = 0*eye(length(Ar)) == Ar*Q_xr + Q_xr*(Ar') + (Br*Q_wr*(Br'));
% 
% [q1, q2, q3, q4] = vpasolve(eq_lyap, [q1, q2, q3, q4])
% 
% Q_xr = [q1 q2;
%         q3 q4];

Q_xr = dlyap(Ar, Br*Q_wr*Br')
Q_yr = Cr*Q_xr*Cr' + Dr*Q_wr*Dr'

sigma_yr = sqrt(Q_yr)
    
%% GLR
% f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])

% False alarm probability:
P_F = 0.0001;

% Probability of missed detection
P_M = 0.01;


h = chi2inv(1 - P_F, 1)/2;                  % Put the threshold from GLR here

% We want to make sure that the threshold P_F = 1 - chi2cdf(2*h, 1) and
% -> h = chi2inv(1 - P_F, 1)/2
% For the window size it's slightly more of a hassle:
% P_D = 
syms zz gg;     % zz is the integrant variable (X), gg is lambda in symbolic

% Density function expression
pd_zz = 1/2*(zz/gg)^(-1/4)*exp((-zz + gg)/2)*besseli(-0.5, sqrt(gg*zz));
p_zz = int(pd_zz, zz,2*h, Inf);  % FILL IN - Integrate over the probability space

eq_1 = p_zz == P_M;  % FILL IN - Equation to be solved
% lambda = double(vpasolve(eq_1, gg)); %1.2625
lambda=1.2625;

%% M
syms M
mu_0 = 0;
mu_1 = f_m(2)*dcgain(RG2(4));

sigma_c = double(Q_yr);

eq_2 = lambda == M*(mu_1 - mu_0)^2/sigma_c^2;

M = double(solve(eq_2, M))

r_window=zeros(floor(M),1);

%%
syms lambda2 M
lambda2 = M*(mu_1 - mu_0)^2/sigma_c^2

pd_zz = 1/2*(zz/lambda2)^(-1/4)*exp((-zz + lambda2)/2)*besseli(-0.5, sqrt(lambda2*zz));
p_zz = int(pd_zz, zz,2*h, Inf);  % FILL IN - Integrate over the probability space

eq_1 = p_zz == P_M;  % FILL IN - Equation to be solved
X = double(vpasolve(eq_1, M, [1, inf])); %1.2625
