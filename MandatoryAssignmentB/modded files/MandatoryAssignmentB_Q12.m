%% 31320 Mandatory Assignment B

MandatoryAssignmentB_Q9

%% Question 12

% Robust stability er frekvens afhængigt så den er defineret som:
% |W_I*T| < 1
% Upper bound of W3:
[magT phase wT] = bode(T, {0.001,1000});

C=squeeze(1./(magT));
sysStab=frd(C,wT);
figure, bodemag(sysStab)
title('$W_3$ upper bound for robust stability','FontName','times','Interpreter','latex')
grid on
saveas(gcf, 'figures/Q12_upper_bound_stability.svg')

W3upperbound=sysStab

save('W3upperbound.mat', 'W3upperbound');


% Robust performance
% Vi har defineret en cirkel til Wp (W1), S og T, men vi mangler WI (W3),
% som kan beregnes ud fra Closed loop og W1
% Slides lecture 16 side 16

% |W1*S| + |W3*T| < 1   -->
% |W3| < (1 - |W1*S|)/|T|

[magW1S,p1]=bode(W1*SenFun,wT);
C=squeeze((1 - magW1S)./(magT));
sysPerf=frd(C,wT);
figure, bodemag(sysPerf)
title('$W_3$ upper bound for robust performance','FontName','times','Interpreter','latex')
grid on
saveas(gcf, 'figures/Q12_upper_bound_performance.svg')


% From upperbound we read that for a static gain W3 needs to be less than
% -23.2 dB

% In other words a high-pass filter could be a good solution.