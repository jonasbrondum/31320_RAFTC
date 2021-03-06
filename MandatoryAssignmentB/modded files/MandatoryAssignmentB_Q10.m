%% 31320 RAFTC
% Question 10

MandatoryAssignmentB_Q9

%% Result validation validate: 
figure ,bode(Pcl('y','r'), (G*K)/(1+G*K),'--r',{0.1,1000})
disp ('validate that G*K/(1+G*K) == transfer function from ref to y')
title('Reference ($r$) to output ($y$)','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$G\cdot K/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on
saveas(gcf, 'figures/Q10_closed_loop_ref2y.svg')

figure ,bode(Pcl('z1','r'), W1/(1+G*K),'--r',{0.1,1000})
disp ('validate that W1/(1+G*K) == transfer function from ref to z1')
title('Reference ($r$) to $z_1$','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$W\_1/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on
saveas(gcf, 'figures/Q10_closed_loop_ref2z1.svg')

figure ,bode(Pcl('z2','r'), (W2*K)/(1+G*K),'--r',{0.1,1000})
disp ('validate that (W2*K)/(1+G) == transfer function from ref to z2')
title('Reference ($r$) to $z_2$','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$W\_2\cdot K/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on
saveas(gcf, 'figures/Q10_closed_loop_ref2z2.svg')

figure ,bode(Pcl('e','r'), 1/(1+G*K),'--r',{0.1,1000})
disp ('validate that 1/(1+G*K) == transfer function from ref to e')
title('Reference ($r$) to $e$','FontName','times','Interpreter','latex')
legend('$P\_{CL}$','$1/(1 + G \cdot K)$','FontName','times','Interpreter','latex')
grid on
saveas(gcf, 'figures/Q10_closed_loop_ref2e.svg')

% Controller assumes positive feedback

hinfnorm(Pcl('y','r'))
hinfnorm(Pcl('z1','r'))
hinfnorm(Pcl('z2','r'))
hinfnorm(Pcl('e','r'))

figure;
step(Pcl('y','r'))
title('Output, $y$, for reference ($r$) step','FontName','times','Interpreter','latex')
saveas(gcf, 'figures/Q10_step_ref2y.svg')



%% Poles and zeros of controller:
[pG, zG] = pzmap(G);

[pK, zK] = pzmap(K);

[pW1, zW1] = pzmap(W1);

[pW2, zW2] = pzmap(W2);

figure;
pzplot(K,G,W1,W2)
title('Pole-zero plot of controller ($K$), plant ($G$), and weight functions $W_1$ and $W_2$','FontName','times','Interpreter','latex')
legend('$K$','$G$','$W\_1$','$W\_2$','Location','southwest','FontName','times','Interpreter','latex')
grid on
saveas(gcf,'figures/Q10_pz_zoomout.svg')

figure;
pzplot(K,G,W1,W2)
title('Pole-zero plot of controller ($K$), plant ($G$), and weight functions $W_1$ and $W_2$','FontName','times','Interpreter','latex')
legend('$K$','$G$','$W\_1$','$W\_2$','Location','southwest','FontName','times','Interpreter','latex')
axis([-5 0.3 -80 80])
grid on
saveas(gcf,'figures/Q10_pz_zoomin.svg')

%% Design verification

figure;
step(c2d(Pcl,T_s))
title('Reference step to y','FontName','times','Interpreter','latex')
saveas(gcf, 'figures/Q10_design_verification.svg')