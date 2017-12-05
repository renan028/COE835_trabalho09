%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular o trabalho 7
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 3     Adaptive parameters
% Caso com observador completo
%----------------------------------------------------------------------

function dx = backstepping_obs(t,x)

global A thetas A0 c1 c2 d1 d2 Gamma gamma kp a w e1 e2 k;

X           = x(1:2); y = e1'*X;
theta       = x(3:5);
lambda      = x(6:7);
eta         = x(8:9);
rho         = x(10);

%% Input
yr = sum(a.*sin(w*t));
dyr = sum(a.*w.*cos(w*t));
ddyr = sum(-a.*w.*w.*sin(w*t));

%% Variables 1
xi = -A0^2 * eta;
Xi = -[A0*eta eta];
v0_1 = lambda(1);
v0_2 = lambda(2);
omega_bar = [0, (Xi(2,:) - y*e1')]';
omega = [v0_2, (Xi(2,:) - y*e1')]';

%% Z
z1 = y - yr;
alpha_bar = -c1*z1 - d1*z1 - xi(2) - omega_bar'*theta;
alpha_1 = rho * alpha_bar;
z2 = v0_2 - rho*dyr - alpha_1;

%% Filtro eta
deta = A0*eta + e2*y;

%% dalpha/dt
dady = rho * (- c1 - d1 + [0,e1']*theta);
dadeta_deta = rho * (e2' * A0^2 * deta + [0,e2'*A0*deta, e2'*eye(2)*deta]*theta);
dadyr = rho*(c1 + d1);
dadtheta = - rho * omega_bar';
dadrho = -(c1 + d1)*z1 - e2'*xi - omega_bar'*theta;

%% Variables 2
tau_1 = (omega - rho*(dyr + alpha_bar)*[e1',0]')*z1;
tau_2 = tau_1 - z2 * (dady * omega);

%% Atualizacao dos parametros
dtheta = Gamma * tau_2;
drho = - gamma * z1 * sign(kp) * (dyr + alpha_bar);
beta = k(2)*v0_1 + dady * (xi(2) + omega'*theta) + ...
    dadeta_deta + dadyr * dyr + (dyr + dadrho) * drho;
u = -c2*z2 + beta + rho*ddyr + dadtheta*dtheta - d2*z2*(dady)^2 - z1*theta(1);

%% Filtros
dlambda = A0*lambda + e2*u;

%% Planta
Phi = [-y 0;0 -y];
F = [e2*u Phi];
dX = A*X + F*thetas;

%% Translation
dx = [dX' dtheta' dlambda' deta' drho]';
