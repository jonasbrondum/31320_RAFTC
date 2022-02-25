clear all;
clc;

% Define the input
syms Vd(t) Vq(t);
ST_input(1,:) = {Vd, Vq};
ST_input(2,:)    = {'V_d', 'V_q'};  % LaTeX % Also as latex(V_d)

% Define the outputs
syms ya(t) yb(t);
ST_meas(1,:)    = {ya, yb};
ST_meas(2,:)    = {'y_1', 'y_2'};   % LaTeX

% Define any parameters
syms R Ld Lq lambda J b P;
ST_par(1,:)    = {R,Ld,Lq,lambda,J,b,P};
ST_par(2,:)    = {'R_s','L_d','L_q','\\lambda_m','J_m','b','P_m'};

% Define the unknown variables
syms ia(t) dia(t) ib(t) dib(t) omega(t) domega(t) T_m(t);
ST_unknowns(1,:) = {ia,dia,ib,dib,omega,domega,T_m};
ST_unknowns(2,:) = {'i_d','\\dot{i}_d','i_q','\\dot{i}_q',...
                    '\\omega','\\dot{\\omega}','T_m'};    % LaTeX

% Define the constraints
ST_cons(1,:) = {'c1','c2','c3','c4','d5','d6','d7','m8','m9'}; % names to be used in the equations
ST_cons(2,:) = {'c_1','c_2','c_3','c_4','d_5','d_6','d_7','m_8','m_9'};     % LaTeX notation  
ST_cons(3,:) = {...
    0 == dia + R/Ld*ia - Lq/Ld*ib*omega - 1/Ld*Vd, ...
    0 == dib + R/Lq*ib + Ld/Lq*ia*omega + lambda/Lq*omega - 1/Lq*Vq, ...
    0 == domega - 1/J*(T_m - b*omega), ...
    0 == 4/(3*P)*T_m - lambda*ib - (Lq - Ld)*ia*ib, ...
    0 == dia - diff(ia,t),...
    0 == dib - diff(ib,t),...
    0 == domega - diff(omega,t),...
    0 == ya - ia,...
    0 == yb - ib};

ST_canfail  = [1 2 3 4 8 9]; % constraints d5, c6 and d7 can't fail
ST_domains  = [1 1 1 1 1 1 1 1 1]; % Number of constraints

% Incidence matrix computed manually
my_ST_IncidenceMatrix = ...
    [ 1 0 0 0  1 1  1 0  1 0 0;
      0 1 0 0  1 0  1 1  1 0 0;
      0 0 0 0  0 0  0 0  1 1 1;
      0 0 0 0  1 0  1 0  0 0 1;
      0 0 0 0 -1 1  0 0  0 0 0
      0 0 0 0  0 0 -1 1  0 0 0
      0 0 0 0  0 0  0 0 -1 1 0
      0 0 1 0  1 0  0 0  0 0 0
      0 0 0 1  0 0  1 0  0 0 0 ];

% Constraints that are non-invertible (write the corresponding "known" variable)
ST_cons_oneway = {[],[],[],[],{ia},{ib},{omega},[],[]};
% Compute incidence matrix
ST_IncidenceMatrix = sa_incidencematrix(ST_cons,ST_input,ST_meas,ST_unknowns,ST_cons_oneway);

% Create the system object
PMSM_sys = sa_create(ST_IncidenceMatrix,ST_cons,ST_input, ST_meas,ST_unknowns, ST_domains,...
           ST_canfail, ST_par);
% Display the system
sa_disp(PMSM_sys);

% Perform a matching
PMSM_sys=sa_match(PMSM_sys,'rank');

% Autogenerate pdf report - requires LaTeX to be installed
sa_report(PMSM_sys,'PMSM_sys','pdf',...
            'analytic_expressions',true);
disp('A report has been generated with the following results:')

% Display the results
disp('Obtained matching:');
sa_disp(PMSM_sys, 't');

disp('Parity relation symbolic form:')
sa_disp(PMSM_sys, 's');
disp('Parity relation analytic form:')
sa_disp(PMSM_sys, 'a');

