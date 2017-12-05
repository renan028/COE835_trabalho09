%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular exemplo 
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 4, 6     Adaptive parameters
%----------------------------------------------------------------------

global Ay By c1 c2 Gamma gamma kp a w;

kp = 1;
Z = [1];
P = [1 2 1];
ss_H = canon(ss(tf(kp*Z,P)), 'companion'); % Planta
[Ay,By,Cy,~,~] = ctrbf(ss_H.A,ss_H.B,ss_H.C);

a = [1 1];
w = [1 1.1];
c1 = 2;
c2 = 2;
Gamma = eye(2);
gamma = 1;
tfinal = 50;

% Initialization
y0  = [5 0]';
theta0 = [0 0]';
p0 = 2;
init = [y0' theta0' p0]';

%% Plots
[T,X] = ode23s('backstepping',tfinal,init,'');

y      = X(:,1);
theta =  X(:,3:4);
p = X(:,5);

yr = a(1)*sin(w(1).*T) + a(2)*sin(w(2).*T);
e =  y - yr;

%Set matlab interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

figure(1)
clf
plot(T,y,T,yr);grid;shg
legend('y','yr','Location','SouthEast')
title('$\epsilon$')
print -depsc2 en3t0
%---------------------------------------------------------------------
