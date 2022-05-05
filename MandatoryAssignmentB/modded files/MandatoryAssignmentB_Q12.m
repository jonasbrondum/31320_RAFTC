%% 31320 Mandatory Assignment B

MandatoryAssignmentB_Q9

%% Question 12

% Maximal multiplicative uncertainty for Robust Stability (RS):
% (BRUGES IKKE) Page 262 about SISO uncertainty




% Robust stability er frekvens afhængigt så den er defineret som:
% |W_I*T| < 1
% Upper bound of W3:
[magT phase wT] = bode(T);

C=squeeze(1./(magT));
sys=frd(C,wT);
figure, bodemag(sys)
title('$W_3$ upper bound','FontName','times','Interpreter','latex')
grid on


% Robust performance
% Vi har defineret en cirkel til Wp (W1), S og T, men vi mangler WI (W3),
% som kan beregnes ud fra Closed loop og W1
% Slides lecture 16 side 16

% |W1*S| + |W3*T| < 1   -->
% |W3| < (1 - |W1*S|)/|T|

[magW1S,p1]=bode(W1*SenFun,wT);
C=squeeze((1 - magW1S)./(magT));
sys=frd(C,wT);
figure, bodemag(sys)
title('Robust performance','FontName','times','Interpreter','latex')
grid on


% From upperbound we read that for a static gain W3 needs to be less than
% -23.2 dB

% In other words a high-pass filter could be a good solution.

%% Test of W3

% High-pass W3
W3 = makeweight(0.05,200,50);

bodemag(W1*SenFun + W3*T)