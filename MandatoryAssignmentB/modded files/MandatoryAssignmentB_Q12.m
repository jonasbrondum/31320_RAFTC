%% 31320 Mandatory Assignment B

MandatoryAssignmentB_Q9

%% Question 12

% Maximal multiplicative uncertainty for Robust Stability (RS):
% (BRUGES IKKE) Page 262 about SISO uncertainty


% Robust stability er vores afstand fra max(CL) til 1

RS_marg = 1 - hinfnorm(CL)


% Robust performance
% Vi har defineret en cirkel til Wp (W1), S og T, men vi mangler WI (W3),
% som kan beregnes ud fra Closed loop og W1
% Slides lecture 16 side 16

% |W1*S| + |W3*T| < 1   -->
% |W3| < (1 - |W1*S|)/|T|

% Upper bound of W3:
W3upperbound = (1 - W1*SenFun)/T;
figure;
sigma(W3upperbound)
title('W3 upper bound')
grid on

% From upperbound we read that for a static gain W3 needs to be less than
% -23.2 dB

% Otherwise a high-pass filter could be a good solution.

%% Test of W3
% Static W3
epsilon = 0.01;
W3 = tf(RS_marg - epsilon,1);
[A,B,C,D] = tf2ss(W3.Numerator{1},W3.Denominator{1});
W3 = ss(A,B,C,D);

% High-pass W3
W3 = makeweight(0.05,200,50);

sigma(W1*SenFun + W3*T)
hinfnorm(W1*SenFun + W3*T)