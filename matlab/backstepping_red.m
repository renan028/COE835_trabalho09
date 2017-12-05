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

function dx = backstepping_red(t,x)

global A thetas N c1 c2 d1 d2 Gamma gamma kp a w e1 e2;

X           = x(1:2); y = e1'*X;
theta       = x(3:5);
lambda      = x(6);
eta         = x(7);
rho         = x(8);

%% Input
yr = sum(a.*sin(w*t));
dyr = sum(a.*w.*cos(w*t));
ddyr = sum(-a.*w.*w.*sin(w*t));

%% Variables 1
xi = -N^2 * eta;
Xi = -[N*eta eta];
v0 = lambda(1);
omega_bar = [0, (Xi - y*e1')]';
omega = [v0, (Xi - y*e1')]';

%% Z
z1 = y - yr;
alpha_bar = -c1*z1 - d1*z1 - xi - omega_bar'*theta + N*y;
alpha_1 = rho * alpha_bar;
z2 = v0 - rho*dyr - alpha_1;

%% Filtro eta
deta = N*eta + y;

%% dalpha/dt
dady = rho * (- c1 - d1 + [0,e1']*theta + N);
dadeta_deta = rho * (N^2 * deta + [0, N*deta, deta]*theta);
dadyr = rho*(c1 + d1);
dadtheta = - rho * omega_bar';
dadrho = -(c1 + d1)*(y - yr) - xi - omega_bar'*theta + N*y;

%% Variables 2
tau_1 = (omega - rho*(dyr + alpha_bar)*[e1',0]')*z1;
tau_2 = tau_1 - z2 * (dady * omega); 

%% Atualiza��o dos par�metros
dtheta = Gamma * tau_2;
drho = - gamma * z1 * sign(kp) * (dyr + alpha_bar);
beta = -N*v0 + dady * (xi + omega'*theta - N*y) + ...
    dadeta_deta + dadyr * dyr + (dyr + dadrho) * drho;
u = -c2*z2 + beta + rho*ddyr + dadtheta*dtheta - d2*z2*(dady)^2 - theta(1)*z1;

%% Filtros
dlambda = N*lambda + u;

%% Planta
Phi = [-y 0;0 -y];
F = [e2*u Phi];
dX = A*X + F*thetas;

%% Translation
dx = [dX' dtheta' dlambda' deta' drho]';    
