%% Mandatory assignment A - 31320 
clear; clc;
%% Question 1:
% 
% Make a structural analysis:
% 1. Determine a complete matching on the unknown variables.
% 2. Find the parity relations in symbolic form.
% 3. Investigate other properties you find relevant from a structural analysis.
% 4. Reformulate the parity relations to an analytic form, as functions only of
% known variables and the system parameters.


% We will use the SA-toolbox to determine most of this

% 1. 


clear;
clc;

% With SA-tool
% Inputs
syms u1(t) u2(t) % Symbolic input declaration
ST_input(1,:) = {u1, u2}; % symbolic variable
ST_input(2,:) = {'u_1\\left(t\\right)','u_2\\left(t\\right)'}; % LaTeX expression


% Measurements
syms y1(t) y2(t) y3(t) % Symbolic measurement declaration
ST_meas(1,:) = {y1, y2, y3}; % symbolic variable
ST_meas(2,:) = {'y_1\\left(t\\right)', 'y_2\\left(t\\right)',...
    'y_3\\left(t\\right)'}; % LaTeX expression


% Parameters
syms J1 J2 J3 k1 k2 b1 b2 b3% Symbolic parameters declaration
ST_parameters(1,:) = {J1 J2 J3 k1 k2 b1 b2 b3}; % symbolic variable
ST_parameters(2,:) = {'J_1', 'J_2', 'J_3', 'k_1', 'k_2,' 'b_1,' 'b_2', 'b_3'}; % LaTeX expression


% Unknowns
syms theta1(t) theta1dot(t) omega1(t) omega1dot(t) theta2(t) theta2dot(t) omega2(t) ...
    omega2dot(t) theta3(t) theta3dot(t) omega3(t) omega3dot(t) d(t) % Symbolic unknowns declaration
ST_unknowns(1,:) = {theta1 theta1dot omega1 omega1dot theta2 theta2dot omega2 omega2dot ...
    theta3 theta3dot omega3 omega3dot d}; % symbolic variable
ST_unknowns(2,:) = {...
	'\\theta_1\\left(t\\right)','\\dot{\\theta_1}\\left(t\\right)',...
    '\\omega_1\\left(t\\right)','\\dot{\\omega_1}\\left(t\\right)',...
	'\\theta_2\\left(t\\right)','\\dot{\\theta_2}\\left(t\\right)',...
    '\\omega_2\\left(t\\right)','\\dot{\\omega_2}\\left(t\\right)',...
    '\\theta_3\\left(t\\right)','\\dot{\\theta_3}\\left(t\\right)',...
    '\\omega_3\\left(t\\right)','\\dot{\\omega_3}\\left(t\\right)',...
    'd\\left(t\\right)'}; % LaTeX expression


% Enter the constraints of the system
% ST_cons(1.:) comprise the Matlab names of constraints
% ST_cons(2.:) comprise the Latex 
% ST_cons(3,:) comprise a cell array list with the expressions of the constraints

ST_cons(1,:) = {'c1','c2','c3','c4','c5','c6','d7','d8','d9','d10','d11','d12','m13','m14','m15'}; % Constraint names
ST_cons(2,:) = {'c_1','c_2','c_3','c_4','c_5','c_6','d_7','d_8','d_9','d_{10}','d_{11}','d_{12}',...
    'm_{13}','m_{14}','m_{15}'}; % Constraint latex names
ST_cons(3,:) = {...
    0 == theta1dot - omega1, ...
    0 == J1*omega1dot - u1 + b1*omega1 + k1*(theta1 - theta2) + d, ...
    0 == theta2dot - omega2, ...
    0 == J2*omega2dot - u2 + b2*omega2 + k1*(theta2 - theta1) + k2*(theta2 - theta3),...
    0 == theta3dot - omega3, ...
    0 == J3*omega3dot + b3*omega3 + k2*(theta3 - theta2),...
    0 == theta1dot - diff(theta1, t),...
    0 == omega1dot - diff(omega1, t),...
    0 == theta2dot - diff(theta2, t),...
    0 == omega2dot - diff(omega2, t),...
    0 == theta3dot - diff(theta3, t),...
    0 == omega3dot - diff(omega3, t),...
    0 == y1 - theta1,...
    0 == y2 - theta2,...
    0 == y3 - theta3};
% NOTE that "diff(.,t)" is used as the differential operator, D == d/dt.


ST_canfail  = [1 2 3 4 5 6 13 14 15]; % constraints d7 - d12 cannot fail
ST_domains  = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
% NOTE NOTE NOTE:
% Incidence matrix MUST have columns ordered as: 
% input, measurements, unknowns
% an x in the book is a "-1" in the software tool
  

% Constraints that are non-invertible (write the corresponding "known" variable)
ST_cons_oneway = {[],[],[],[],[],[],{theta1},{omega1},{theta2},{omega2},{theta3},{omega3},...
    [],[],[]};
% Automatic incMatrix generation
ST_IncMat = sa_incidencematrix(ST_cons,...
                                    ST_input,ST_meas,...
                                    ST_unknowns,ST_cons_oneway);
						
								
ST_sys =...
    sa_create(ST_IncMat,ST_cons,...
    ST_input, ST_meas,ST_unknowns,...
    ST_domains, ST_canfail,...
	ST_parameters);



sa_disp(ST_sys);

% Perform a matching
ST_sys=sa_match(ST_sys,'rank');

% Autogenerate pdf report - requires LaTeX to be installed
sa_report(ST_sys,'MandatoryAssignmentA_Q1','pdf',...
            'analytic_expressions',true);
disp('A report has been generated with the following results:')

% Display the results
disp('Obtained matching:');
sa_disp(ST_sys, 't');

disp('Parity relation symbolic form:')
sa_disp(ST_sys, 's');
disp('Parity relation analytic form:')
sa_disp(ST_sys, 'a');
