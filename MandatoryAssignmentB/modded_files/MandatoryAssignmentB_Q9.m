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
system = tf(num(3,:),den)

bode(system)

[K, CL, gamma] = hinfsyn(system, 1, 1)

%% Test of implementation:
% Augmented plant for synthesis
%                                       +------+
%            +--------------------------|W_p(s)|----> z1
%            |                          +------+
%            |                        +------+
%            |                 +------|W_u(s)|------> z2
%            |                 |      +------+
%         +  | e +------+  u(s)| +------+             
% ref ---->O-----| K(s) |--------| G(s) |-----------> y
%         -|     +------+        +------+       |         
%          |                                    |
%          +------------------------------------+
%                  
%
% UNFINISHED:
[A,B,C,D] = tf2ss(k,[Is 0 0]);
G = ss(A,B,C,D);
Ks = Kc*((s+tau1)/s)*((tau2*s+1)/(tau3*s+1));

M = 2;
A = 0.01;
wb = pi/5;
W_p = (s/M + wb)/(s + wb*A) ;% Choose "control sensitivity weight" as highpass filter
[Awp,Bwp,Cwp,Dwp] = ssdata(W_p);
Wp = ss(Awp,Bwp,Cwp,Dwp);

%ny = 2;
%nu = 2;

G.InputName = 'ug';
G.OutputName= 'yg';

Wp.InputName = 'uWp';
Wp.OutputName = 'yWp';

Ks.InputName = 'uKs';
Ks.OutputName = 'yKs';

blksys = append(G,Wp,Ks);

% specify connections
S1 = sumblk('ug = yKs',1);
S2 = sumblk('y = yg',1);
S3 = sumblk('uWp = e',1);
S4 = sumblk('uKs= e',1);
S5 = sumblk('zp = yWp',1); 
S6 = sumblk('e = r - yg',1); 
S = append(S1,S2,S3,S4,S5,S6);
 


 % connect to obtain generalised plant
% P    = [Wp , Wp*G, Wp*G;
%         0  , 0   , Wu;
%         I  , G   , G];
Pcl  = connect(blksys,S,{'r'},{'y','zp','e'});

Gcl = feedback(Ks*G,1);
% validate: 
figure ,bode(Pcl('y','r'), (G*Ks)/(1+G*Ks),'--r')
disp ('validate that G*Ks/(1+G*Ks) == transfer function from ref to y')
figure ,bode(Pcl('zp','r'), Wp/(1+G*Ks),'--r')
disp ('validate that Wp/(1+G*Ks) == transfer function from ref to zp')
figure ,bode(Pcl('e','r'), 1/(1+G*Ks),'--r')
disp ('validate that 1/(1+G*Ks) == transfer function from ref to e')

% Controller assumes positive feedback

hinfnorm(Pcl('y','r'))
hinfnorm(Pcl('zp','r'))
hinfnorm(Pcl('e','r'))