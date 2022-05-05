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

s = tf ('s')

%notch filter

%According to Brian Douglas, i want to lower the peak with 0.001


%The two notches match the eigenvalues of G
wb=31.67;
Q=1;

W0 = (s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)

wb=65.8;
Q=1;

W0 = W0*(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)


[A,B,C,D] = tf2ss(W0.Numerator{1},W0.Denominator{1});

notches=ss(A,B,C,D)

W1 = notches*makeweight(2,10,0.01);

W1=W3+1

%W1=1/feedback(1,G)*2





W2= makeweight(0.1,8,10)%setup 1

W2= makeweight(0.8,8,100)%setup 3

W2= makeweight(0.8,25,100)%setup 5

%W2=notches

[K,CL,gamma] = mixsyn(G,W1,W2,W3); %Last argument makes the function try to force gamma to 1
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

figure
sigma(feedback(1,G),feedback(1,G)*W1)



%% Simulink 

[numK, denK]=ss2tf(K.A, K.B, K.C, K.D,1);

%[numW1, denW1]=ss2tf(W1.A, W1.B, W1.C, W1.D,1);


%[numW2, denW2]=ss2tf(W2.A, W2.B, W2.C, W2.D,1);


%Further, it is important that the controller does not have too high gains.
%Our controller has crazy high gains


x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.0004;                    % Sampling period


