clear all;
close all;
clc;

%% load real data

%load('Experiment\first_run.mat');
load('Experiment\Second_run.mat');



u_ex = timetable([zeros(size(u(:,2))) u(:,2)],'SampleRate', 1/T_s);

%u_ex = timetable([zeros(size(u_unsat(:,2))) u_unsat(:,2)],'SampleRate', 1/T_s);
%u_ex = timetable([zeros(size(u_unsat(:,2))) u_unsat(:,2)],'SampleRate', 1/T_s);

y_ex = timetable(y(:,2:4),'SampleRate', 1/T_s);

%% Filtering

s = tf('s')

sys = -s^2-(b_2*s-k_2-k_1)/J_2
zeta=0.7;
omega=sqrt(k_2/J_3)
% Originally we chose omega = 10/4
omega=10; % From experiment 28/02 -> omega = 20 was a good value

%omega=1; % From experiment 28/02 -> omega = 20 was a good value

omega=20; % From experiment 28/02 -> omega = 20 was a good value

lowpass1=omega^2/(s^2+2*zeta*omega*s+omega^2)


%% Discretize

s = tf('s')


F=[ 0 1/J_2 k_1/J_2 (-s^2*J_2-b_2*s-k_1-k_2)/J_2 k_2/J_2;
    0 0 0 k_2/J_3 (-J_3*s^2-b_3*s-k_2)/J_3];
    
%Q = [lowpass1 0;
%    0 lowpass1]

RG1 = c2d(lowpass1*F(1,:),T_s,'tustin')
RG2 = c2d(lowpass1*F(2,:),T_s,'tustin')

[num1 den1] = tfdata(RG1);
[num2 den2] = tfdata(RG2);


num12=cell2mat(num1(2));
num13=cell2mat(num1(3));
num14=cell2mat(num1(4));
num15=cell2mat(num1(5));


den12=cell2mat(den1(2));
den13=cell2mat(den1(3));
den14=cell2mat(den1(4));
den15=cell2mat(den1(5));



num24=cell2mat(num2(4));
num25=cell2mat(num2(5));
den24=cell2mat(den2(4));
den25=cell2mat(den2(5));





simTime = 45;                   % Simulation duration in seconds
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
f_m = [0,0,0];
u_fault = 0;                    % Disable VA meachanism
detect_time = 45;
f_u_time=45;
%% Plot experiment data

load('Experiment\first_run.mat');
r3=(J_2*k_2)/(J_3*k_1+J_3*k_2)*r(:,2)+r(:,3);

figure
hold on
plot(r(:,1),r(:,2),'LineWidth',2)
plot(r(:,1),r(:,3),'LineWidth',2)
plot(r(:,1),r3,'LineWidth',2)
legend('RG1', 'RG2', 'RG3');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{r}(t)$','FontName','times','FontSize',16,'Interpreter','latex')
hold off



load('Experiment\Second_run.mat');
r3=(J_2*k_2)/(J_3*k_1+J_3*k_2)*r(:,2)+r(:,3);

figure
hold on
plot(r(:,1),r3,'LineWidth',2)
plot(r(:,1),r(:,2),'LineWidth',2)
plot(r(:,1),r(:,3),'LineWidth',2)
legend('RG3', 'RG1', 'RG2');
xlabel('Time [sec]','FontName','times','FontSize',16,'Interpreter','latex')
ylabel('$\mathbf{r}(t)$','FontName','times','FontSize',16,'Interpreter','latex')
set(l,'FontName','times','FontSize',16,'Interpreter','latex');
hold off


%Residual 3

