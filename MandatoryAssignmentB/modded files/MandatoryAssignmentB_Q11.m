%% 31320 Mandatory Assignment B
% Question 11
MandatoryAssignmentB_Q9


sim('threeDiskOscillatorRig_noFric');

%% Plotting position

figure
hold on
plot(position,'LineWidth',2)
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$r(t)$','FontName','times','FontSize',16,'Interpreter','latex')
saveas(gcf,'figures/Q11_position.svg')


%% Plotting control signal
figure
hold on
plot(input,'LineWidth',2)
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$u(t)$','FontName','times','FontSize',16,'Interpreter','latex')
saveas(gcf,'figures/Q11_input.svg')

