%% 31320 RAFTC
% Question 10

MandatoryAssignmentB_Q9

%% Result validation validate: 
figure ,bode(Pcl('y','r'), (G*K)/(1+G*K),'--r',{0.1,1000})
disp ('validate that G*K/(1+G*K) == transfer function from ref to y')
title('Reference ($r$) to output ($y$)','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$G\cdot K/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on

figure ,bode(Pcl('z1','r'), W1/(1+G*K),'--r',{0.1,1000})
disp ('validate that W1/(1+G*K) == transfer function from ref to z1')
title('Reference ($r$) to $z_1$','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$W_1/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on

figure ,bode(Pcl('z2','r'), (W2*K)/(1+G*K),'--r',{0.1,1000})
disp ('validate that (W2*K)/(1+G) == transfer function from ref to z2')
title('Reference ($r$) to $z_2$','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$W_2\cdot K/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on

figure ,bode(Pcl('e','r'), 1/(1+G*K),'--r',{0.1,1000})
disp ('validate that 1/(1+G*K) == transfer function from ref to e')
title('Reference ($r$) to $e$','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$1/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on

% Controller assumes positive feedback

hinfnorm(Pcl('y','r'))
hinfnorm(Pcl('z1','r'))
hinfnorm(Pcl('z2','r'))
hinfnorm(Pcl('e','r'))

figure;
step(Pcl('y','r'))
title('Output, $y$, for reference ($r$) step','FontName','times','Interpreter','latex')



%% Poles and zeros of controller:
pzG = pzmap(G);

pzK = pzmap(K);

pzW1 = pzmap(W1);

pzW2 = pzmap(W2);
