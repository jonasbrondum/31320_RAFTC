%MandatoryAssignmentB_Q13

W3=tf([0.833 0],[1 0.089]);


%The given weight function for W_I is just above the lower bound between
%4.8e-07 rad/s and 23.6 rad/s. It is a good approximation!


figure
hold on
title('W_I and the lower bound function');
sigma(W3)
sigma(W3lowerbound)
hold off


%% Redefinition of weights for Simulink simulation

% First we need to make some weights relevant as to the design
% requirements:

% We want gamma ~= 1
% And we don't want too high gains
% We want (almost) integral action meaning (almost) no steady state error
% We don't want the actuator to go into saturation for a step-response

% Example from book
% Uses the Robust Control toolbox



%It is very hard to get a Robust system that is fast. Either we go slow for
%stability, or we go fast and get oscilations.

%Fast setup, not very robust stable:

%Setup 1
%M=5; wb=2; A=1.e-4; % Hvad betyder disse og hvor i bogen kommer de fra?


%M=5; wb=2; A=1.e-4;%Setup 2



M=60; wb=2; A=1.e-5;%Setup 1
W1 = tf([1/M wb], [1 wb*A]);


[A,B,C,D] = tf2ss(W1.Numerator{1},W1.Denominator{1});

W1 = ss(A,B,C,D);

%Perhaps W1 should be a notch?


wb=0.1;%Setup 1
zeta=10

W1=tf([0 2*zeta*wb 0], [1 2*zeta*wb wb^2])

%W1=tf([1 0 wb^2], [1 2*zeta wb^2])
W1=1-W1

[A,B,C,D] = tf2ss(W1.Numerator{1},W1.Denominator{1});

W1 = ss(A,B,C,D);



%W1= makeweight(50,0.05,0.8)%Setup 2

%W2= makeweight(0.1,0.1,20)%Old, semi-stable "fast", Setup 2






%Because S + T = I, mixsyn cannot make both S and T small 
%(less than 0 dB) in the same frequency range. 
%Therefore, when you specify weights for loop shaping, 
%there must be a frequency band in which both W1 and W3 are below 0 dB.

%W1= makeweight(50,0.05,0.8)%Robust stable, Setup 3

%W2= makeweight(0.1,0.1,20)%Robust stable, Setup 3


%Setup 2 gives us slight osscilations at 4 Hz, or 25 radians per second.
%We know our approximation of W3 becomes unstable outside 25 rad/s, so
%perhaps this is as good as it gets?
%W2= makeweight(0.01,5,20)%Setup 2
W2= makeweight(0.1,8,10)%setup 1




[K,CL,gamma] = mixsyn(G,W1,W2,W3, 1); %Last argument makes the function try to force gamma to 1
% [K,CL,gamma] = mixsyn(G,W1,[],[]);
gamma

% First, compare the resulting sensitivity S and complementary sensitivity 
% T to the corresponding weighting functions W1 and W3. 

L = G*K;
I = eye(size(L));
SenFun = feedback(I,L); 
T= I-SenFun;

close all;

figure;
sigma(SenFun,'b',W1,'b--',T,'r',W2,'r--',{0.01,1000})
legend('S','W1','T','W3')

figure;
sigma(L,'b',W1,'r--',1/W2,'g--',{0.01,1000})
legend('L','W1','1/W2')

%Singular values look good, gamma is very close to 1.

sigma(SenFun,'b',L,'r',T,'g',gamma/W1,'b-.',ss(gamma/W2),'r-.',gamma/W3,'g-.',{1e-3,1e3})
legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest')
grid

%Check phasemargin
%Marg = allmargin(G*K)

%%




RS_Marg= 1-hinfnorm(CL)

%Our new controller has a higher robust stability margin, since the closed
%loop is significantly larger than in Q12

%The robust performance
hinfnorm(W1*SenFun+W3*T)



figure
hold on
title('W_I and the lower bound function');
sigma(W1*SenFun+W3*T)
hold off


figure
hold on
title('W_I and the lower bound function');
sigma(W3)
sigma(W3lowerbound)

hold off

figure
title('W1, W2, W3')
bodemag(W1,W2,W3)


%% Simulink 

[numK, denK]=ss2tf(K.A, K.B, K.C, K.D,1);

[numW1, denW1]=ss2tf(W1.A, W1.B, W1.C, W1.D,1);


[numW2, denW2]=ss2tf(W2.A, W2.B, W2.C, W2.D,1);



x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.0004;                    % Sampling period


