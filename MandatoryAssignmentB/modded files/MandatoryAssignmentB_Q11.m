%% 31320 Mandatory Assignment B
% Question 11
MandatoryAssignmentB_Q9


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

