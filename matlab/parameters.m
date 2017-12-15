clear;
clc;
close all;
global e1 e2;

% PRINT = true;
PRINT = false;

%Simulation time
tfinal = 100;

% Unit vectors
e1 = [1 0]';
e2 = [0 1]';

%Obsever
s = tf('s');
A0 = s + 1;

%% First parameters

%Plant
kp_1 = 2;
a1_1 = 2;
a0_1 = 1;

%Model
km_1  = 1;
am1_1 = 4;
am0_1 = 4;

%Reference
a_1 = [1 1];
w_1 = [1 3];

%Initial conditions
X0_1  = [0 0]';
X0m_1  = [0 0]';
w10_1 = 0;
w20_1 = 0;
Psi0_1 = [0;0;0;0];
xi0_1 = 0;
Omega0_1 = [0;0;0;0];
rho0_1 = 1;

%Adaptation gain and filter ctes
N_1 = -1;
Gamma_1 = 1;
gamma_1 = 1;
c1_1 = 1;
c2_1 = 1;
d1_1 = 1;
d2_1 = 1;

%% Second parameters

%Plant
kp_2 = 2;
a1_2 = 2;
a0_2 = 1;

%Model
km_2  = 1;
am1_2 = 4;
am0_2 = 4;

%Reference
a_2 = [1 1];
w_2 = [1 3];

%Initial conditions
X0_2  = [0 0]';
X0m_2  = [0 0]';
w10_2 = 0;
w20_2 = 0;
Psi0_2 = [0;0;0;0];
xi0_2 = 0;
Omega0_2 = [0;0;0;0];
rho0_2 = 1;

%Adaptation gain and filter ctes
N_2 = -1;
Gamma_2 = 1;
gamma_2 = 1;
c1_2 = 1;
c2_2 = 1;
d1_2 = 1;
d2_2 = 1;
