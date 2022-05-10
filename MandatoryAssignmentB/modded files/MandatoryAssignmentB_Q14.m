MandatoryAssignmentB_Q13

W3=tf([0.833 0],[1 0.089]);


%The given weight function for W_I is just above the lower bound between
%4.8e-07 rad/s and 23.6 rad/s. It is a good approximation!
%%

figure
hold on
title('W_I and the lower bound function');
%sigma(W3upperbound, {0.001,10^3})
sigma(W3lowerbound, {0.001,10^3})
sigma(W3,{0.001,10^3})
%ylim([-30,50]);
grid on
title('Lower boundary and the given $W_3$ transfer function','FontSize',16,'Interpreter','latex');
legend('$W\_I$','Lower boundary','FontName', 'times','Interpreter','latex')

saveas(gcf,'figures\Q14_lowerbound_WI.svg')

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
M=10; wb=3; A=1.e-4; % Hvad betyder disse og hvor i bogen kommer de fra?

%M=5; wb=5; A=1.e-3; % Hvad betyder disse og hvor i bogen kommer de fra?

%1/M is the gain after the cutoff. We need a lower gain in the second
%iteration, as we have a lower

W1 = tf([1/M wb], [1 wb*A]);


[A,B,C,D] = tf2ss(W1.Numerator{1},W1.Denominator{1});

W1 = ss(A,B,C,D);

W2= makeweight(0.2,40,100)%setup 1




%wb=25;


%W0 = W0*(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)





%Because S + T = I, mixsyn cannot make both S and T small 
%(less than 0 dB) in the same frequency range. 
%Therefore, when you specify weights for loop shaping, 
%there must be a frequency band in which both W1 and W3 are below 0 dB.

%W1= makeweight(50,0.05,0.8)%Robust stable, Setup 3

%W2= makeweight(0.1,0.1,20)%Robust stable, Setup 3


%Setup 2 gives us slight osscilations at 4 Hz, or 25 radians per second.
%We know our approximation of W3 becomes unstable outside 25 rad/s, so
%perhaps this is as good as it gets?





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

figure
hold on
title('W1, W2, W3')
bodemag(W1,W2,W3, {0.01, 1000})
grid on
title('Weight functions for $G_p$','FontSize',16,'Interpreter','latex')
legend('$W\_p$','$W\_u$','$W\_I$','FontName', 'times','Interpreter','latex')
saveas(gcf,'figures,Q14_weights.svg')
hold off


%% Stability conditions

figure
hold on
title('Nominal performance');
bodemag(W1*SenFun)
hold off


figure
hold on
title('Robust stability');
bodemag(W3*T)
hold off



figure
hold on
bodemag(W1*SenFun)
bodemag(W3*T)
bodemag(W1*SenFun+W3*T)
legend('$|W\_p \cdot S|$','$|W\_I \cdot T|$','$|W\_p \cdot S| + |W\_I \cdot T|$','FontName', 'times','Interpreter','latex')
%legend('W1*S', 'W3*T', 'W1*S+W3*T')
title('Multiplicative uncertainty conditions');
grid on
ylim([-30,5])
xlim([0.01, 10^2])
saveas(gcf,'figures\Q14_perfomance.svg')
hold off



%% Simulink 

[numK, denK]=ss2tf(K.A, K.B, K.C, K.D,1);

[numW1, denW1]=ss2tf(W1.A, W1.B, W1.C, W1.D,1);


[numW2, denW2]=ss2tf(W2.A, W2.B, W2.C, W2.D,1);


%Further, it is important that the controller does not have too high gains.
%Our controller has crazy high gains


x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.0004;                    % Sampling period

%%
L1 = G*1.13*K; % Ziegler-Nichols Controller
T1 = feedback(L1,1);
[Smarg1,Dstab1,Report1] = robuststab(T1)

[Smarg1,Dstab1,Report1] = robuststab(G)
