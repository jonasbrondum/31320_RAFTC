%% 31320 Mandatory Assignment B
% Question 11
MandatoryAssignmentB_Q9

%% Sampling period

T_s = 0.004;                    % Sampling period
%% Redefinition of weights for Simulink simulation

% First we need to make some weights relevant as to the design
% requirements:

% We want gamma ~= 1
% And we don't want to high gains
% We want (almost) integral action meaning (almost) no steady state error
% We don't want the actuator to go into saturation for a step-response

% Example from book
% Uses the Robust Control toolbox
M=0.2; wb=32; A=1.e-4;
W1 = tf([1/M wb], [1 wb*A]);

[A,B,C,D] = tf2ss(W1.Numerator{1},W1.Denominator{1});

W1 = ss(A,B,C,D);


[K,CL,gamma] = mixsyn(G,W1,W2ss,[], 20); %Last argument makes the function try to force gamma to 1
% [K,CL,gamma] = mixsyn(G,W1,[],[]);
gamma

% First, compare the resulting sensitivity S and complementary sensitivity 
% T to the corresponding weighting functions W1 and W3. 

L = G*K;
I = eye(size(L));
S = feedback(I,L); 
T= I-S;

close all;

figure;
sigma(S,'b',W1,'b--',T,'r',W2,'r--',{0.1,1000})
legend('S','W1','T','W3')

figure;
sigma(L,'b',W1,'r--',1/W2,'g--',{0.1,1000})
legend('L','W1','1/W2')

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
%
% UNFINISHED:
% [A,B,C,D] = tf2ss(k,[Is 0 0]);
% G = ss(A,B,C,D);
% Ks = Kc*((s+tau1)/s)*((tau2*s+1)/(tau3*s+1));
% 
% M = 2;
% A = 0.01;
% wb = pi/5;
% W_p = (s/M + wb)/(s + wb*A) ;% Choose "control sensitivity weight" as highpass filter
% [Awp,Bwp,Cwp,Dwp] = ssdata(W_p);
% Wp = ss(Awp,Bwp,Cwp,Dwp);

%ny = 2;
%nu = 2;

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
% validate: 
figure ,bode(Pcl('y','r'), (G*K)/(1+G*K),'--r',{0.1,1000})
disp ('validate that G*K/(1+G*K) == transfer function from ref to y')
figure ,bode(Pcl('z1','r'), W1/(1+G*K),'--r',{0.1,1000})
disp ('validate that W1/(1+G*K) == transfer function from ref to z1')
figure ,bode(Pcl('z2','r'), W2/(1+G*K),'--r',{0.1,1000})
disp ('validate that W2/(1+G*K) == transfer function from ref to z2')
figure ,bode(Pcl('e','r'), 1/(1+G*K),'--r',{0.1,1000})
disp ('validate that 1/(1+G*K) == transfer function from ref to e')

% Controller assumes positive feedback

hinfnorm(Pcl('y','r'))
hinfnorm(Pcl('z1','r'))
hinfnorm(Pcl('z2','r'))
hinfnorm(Pcl('e','r'))

figure;
step((Pcl))
figure;

step(c2d(Pcl,T_s))


%% Simulation in Simulink
%Simulink outputs very different results from Henriks method. This is
%probably due to his system being continous, while the simulink model is
%discrete


[numK, denK]=ss2tf(K.A, K.B, K.C, K.D,1)

[numW1, denW1]=ss2tf(W1.A, W1.B, W1.C, W1.D,1)


[numW2, denW2]=ss2tf(W2ss.A, W2ss.B, W2ss.C, W2ss.D,1)



x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.0004;                    % Sampling period
%[numW1 denW1] = (tfdata(W1));

%numW1=cell2mat(numW1)
%denW1=cell2mat(denW1)


%W2 = tf(0.1,1);

%[numW2 denW2] = cell2mat(tfdata(W2));

%[numKtf denKtf] = tfdata(Ktf);




%% Example from book
% Uses the Robust Control toolbox
%G=tf(200,conv([10 1],conv([0.05 1],[0.05 1]))); % Plant is G.
M=1.5; wb=10; A=1.e-4;
Wp = tf([1/M wb], [1 wb*A]); Wu = 1; % Weights.
% Find H-infinity optimal controller:
[khinf,ghinf,gopt] = mixsyn(G,Wp,Wu,[]);
Marg = allmargin(G*khinf) % Gain and phase margins



