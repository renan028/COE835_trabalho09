%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular o trabalho 8
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 3     Adaptive parameters
% Caso com observador de ordem reduzida
%----------------------------------------------------------------------

function dx = backstepping_red_direct(t,x)

global km am1 am0 N c1 c2 d1 d2 Gamma gamma kp a1 a0 a w Af Bf e1;

X           = x(1:2); y = e1'*X;
Xm          = x(3:4); ym = e1'*Xm;
w1          = x(5);
w2          = x(6);
Psi         = x(7:10);
xi          = x(11);
Omega       = x(12:15);
rho         = x(16);

%% Input
r = sum(a.*sin(w*t));

%% Output error
z1 = y - ym;

%% Filters
dxi = N*xi -(am0 + am1*N + N^2)*z1 - km*r;
v0 = Omega(1);
Omega_bar = Omega; Omega_bar(1) = 0;
dOmega_bar = N*Omega_bar + km*[0; -[w1; y; w2]];

%% Stabilization Function
alpha_bar = (-c1 + d1/N + am1 + N)*z1 - xi - Omega_bar'*Psi;
alpha = rho * alpha_bar;

% Partial derivatives
dadz1 = rho * (- c1 + d1/N + am1 + N);
dadrho = alpha_bar;
dadxi = -rho;
dadPsi = - rho * Omega_bar';
dadOmega_bar = - rho * Psi';

%% z2
z2 = v0 - alpha;

%% rho update law
drho = -gamma*sign(kp/km)*alpha_bar*z1;

%% Axuiliary vars
beta = -N*v0 + dadz1 * (xi + Omega'*Psi - (am1 + N)*z1) + dadrho*drho + dadxi*dxi + dadOmega_bar*dOmega_bar;
% beta = -N*v0 + dadz1 * (xi + Omega'*Psi - (am1 + N)*z1) + dadrho*drho + dadxi*dxi;

%% Tuning functions
tau_1 = (Omega - rho*alpha_bar*[e1;0;0])*z1;
tau_2 = tau_1 - z2 * (dadz1 * Omega);

%% Parameter Update
dPsi = Gamma * tau_2;
u = 1/km*(-c2*z2 - Psi(1)*z1 + beta + dadPsi*dPsi + d2/N*z2*dadz1^2);

%% Derivatives
dX = [-a1 1; -a0 0]*X + [0;kp]*u;
dXm = [-am1 1;-am0 0]*Xm + [0;km]*r;
dv0 = N*v0 + km*u;
dOmega = [dv0 dOmega_bar(2:end)']';
dw1 = Af*w1 + Bf*u;
dw2 = Af*w2 + Bf*y;

dx = [dX' dXm' dw1' dw2' dPsi' dxi' dOmega' drho']';
