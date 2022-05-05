
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

