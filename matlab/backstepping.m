%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular exemplo 
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 4, 6     Adaptive parameters
% Caso sem observador
%----------------------------------------------------------------------

function dx = backstepping(t,x)

global Ay By c1 c2 Gamma gamma kp a w;

y      = x(1:2);
theta  = x(3:4);
p      = x(5);

%% Input
yr = a(1) * sin(w(1)*t) + a(2) * sin(w(2)*t);
dyr = a(1) * w(1) * cos(w(1)*t) + a(2) * w(2) * cos(w(2)*t);
ddyr = -a(1) * w(1)^2 * sin(w(1)*t) - a(2) * w(2)^2 * sin(w(2)*t);

%% Variáveis de controle
phi_1 = [0 -y(1)]';
phi_2 = [-y(1) 0]';
%% Z
alpha_1 = -c1 * (y(1) - yr) - phi_1'*theta; 
z1 = y(1) - yr;
z2 = y(2) - alpha_1 - dyr; 
tau_1 = phi_1 * z1;
tau_2 = tau_1 + (phi_2 + c1 * phi_1) * z2; 

%% Atualização dos parâmetros
dtheta = Gamma * tau_2;
alpha_2 = -c2 * z2 - z1 - theta'*(phi_2 + c1*phi_1) - c1*y(2) - ...
    phi_1'*dtheta + c1*dyr;
u_bar = alpha_2 + ddyr;
u = p * u_bar;
dp = - gamma * sign(kp) * u_bar * z2;


%% Planta
dy = Ay*y + By*u;

%% Translation
dx = [dy' dtheta' dp]';    
