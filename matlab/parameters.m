clear;
clc;
close all;
global A e1 e2;

PRINT = true;
%PRINT = false;

%Simulation time
tfinal = 100;

% Unit vectors
e1 = [1 0]';
e2 = [0 1]';

% System matrix
A = [0 1;0 0];

%% First parameters

kp_1 = 5;
Z_1 = [1];
P_1 = [1 2 1];
thetas_1 = [kp_1 P_1(2) P_1(3)]';

k_1 = [1 1]';

%Initial conditions
X0_1  = [0 0]';
theta0_1 = [0 0 0]';
lambda0_1 = [0 0]';
eta0_1 = [0 0]';
rho0_1 = 1;

%Adaptation gain
Gamma_1 = 1;
gamma_1 = 1;
c1_1 = 1;
c2_1 = 1;
d1_1 = 1;
d2_1 = 1;

%Reference
a_1 = [1 1];
w_1 = [1 3];

%% Second parameters

kp_2 = 2;
Z_2 = [1];
P_2 = [1 -2 1];
thetas_2 = [kp_2 P_2(2) P_2(3)]';

k_2 = [1 1]';

%Initial conditions
X0_2  = [10 0]';
theta0_2 = [0 0 0]';
lambda0_2 = [0 0]';
eta0_2 = [0 0]';
rho0_2 = 1;

%Adaptation gain
Gamma_2 = 10;
gamma_2 = 1;
c1_2 = 1;
c2_2 = 1;
d1_2 = 1;
d2_2 = 1;

%Reference
a_2 = [1 2];
w_2 = [1 5];
