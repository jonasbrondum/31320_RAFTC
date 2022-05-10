
MandatoryAssignmentB_Q14
J_1 = 0.0325;  %CHANGED

sim('threeDiskOscillatorRig_noFric');


%% Plotting with nominal controller

load('nominal');

figure
hold on
plot(nominalPosition,'LineWidth',2)
plot(position,'LineWidth',2)
title('Step response of the two systems');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('Disc 3 position [rad]','FontName','times','FontSize',16,'Interpreter','latex')
legend('$G$ Nominal system','$G_p$ Uncertain system','FontName', 'times','Interpreter','latex')
saveas(gcf,'figures\Q15_position.svg')
hold off

figure
hold on
plot(nominalInput,'LineWidth',2)
plot(input,'LineWidth',2)
title('Control signal of the two systems');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('Input signal [V]','FontName','times','FontSize',16,'Interpreter','latex')
legend('$G$ Nominal system','$G_p$ Uncertain system','FontName', 'times','Interpreter','latex')

saveas(gcf,'figures\Q15_input.svg')
hold off


%%



%Nominal system

load('ECP_values.mat');
% Physical system parameters
J_1 = ECP_values(1);            % Disk 1 inertia kgm^2

sim('threeDiskOscillatorRig_noFric');

%% PLotting position

figure
hold on
plot(position,'LineWidth',2)
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{r}(t)$','FontName','times','FontSize',16,'Interpreter','latex')
hold off

%% Plotting control signal
figure
hold on
plot(input,'LineWidth',2)
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{r}(t)$','FontName','times','FontSize',16,'Interpreter','latex')
hold off


%% Getting the coloumb friction bandwitdh
% 
% mdl = 'coloumb_bandwitdh';
% load_system(mdl);
% io(1) = linio('coloumb_bandwitdh/Gain1',1,'input');
% io(2) = linio('coloumb_bandwitdh/Gain',1,'output');
% linsys1 = linearize(mdl,io);
% sys = tf(linsys1)
% bode(sys)



