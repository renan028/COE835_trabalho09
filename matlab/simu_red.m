%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular o trabalho 8 
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 3     Adaptive parameters
% Com observador de ordem reduzida
%----------------------------------------------------------------------

global km am1 am0 N c1 c2 d1 d2 Gamma gamma kp a1 a0 dc a w Af Bf;

sim_str = strcat('');

%% ------------------------------------------------- Simulation 1 (default)

%Plant
kp = kp_1;
a1 = a1_1;
a0 = a0_1;
P = tf(kp,[1 a1 a0]);

%Model
km = km_1;
am1 = am1_1;
am0 = am0_1;
Pm = tf(km,[1 am1 am0]);

%2DOF Control ideal parameters and lambda filter
[t1, tn, t2, t2n, L] = find2DOFparameters(P,Pm,A0);
Psis_1 = 1/t2n*[1 t1 tn t2];

%u and y filter
ss_f = canon(ss(tf(1,L)), 'companion');
Af = ss_f.A';
Bf = ss_f.C';

% Reference
dc = dc_1;
a = a_1;
w = w_1;

%Initialization
X0  = X0_1;
X0m  = X0m_1;
w10 = w10_1;
w20 = w20_1;
Psi0 = Psi0_1;
xi0 = xi0_1;
Omega0 = Omega0_1;
rho0 = rho0_1;

init = [X0' X0m' w10' w20' Psi0' xi0' Omega0' rho0]';

%Adaptation gain and filter ctes
N = N_1;
c1 = c1_1;
c2 = c2_1;
d1 = d1_1;
d2 = d2_1;
Gamma = Gamma_1*eye(length(Psi0));
gamma = gamma_1;

%Simulation
[T_1,X_1] = ode23s('backstepping_red_direct',tfinal,init,'');
y_1       = X_1(:,1);
ym_1      = X_1(:,3);
Psi_1     = X_1(:,7:10);
rho_1     = X_1(:,end);
tilPsi_1  = -Psi_1 + Psis_1.*ones(length(Psi_1),length(Psis_1));
modPsi_1  = sqrt(sum(Psi_1.^2,2));
e0_1 =  y_1 - ym_1;

r_1 = 0;
for i=1:length(a)
    r_1 = r_1 + a(i)*sin(w(i).*T_1);
end

%% --------------------------------------------------- Simulation 2 (Gamma)
changed = 1;

%Plant
kp = kp_1;
a1 = a1_1;
a0 = a0_1;
P = tf(kp,[1 a1 a0]);

%Model
km = km_1;
am1 = am1_1;
am0 = am0_1;
Pm = tf(km,[1 am1 am0]);

%2DOF Control ideal parameters and lambda filter
[t1, tn, t2, t2n, L] = find2DOFparameters(P,Pm,A0);
Psis_2 = 1/t2n*[1 t1 tn t2];

%u and y filter
ss_f = canon(ss(tf(1,L)), 'companion');
Af = ss_f.A';
Bf = ss_f.C';

% Reference
dc = dc_1;
a = a_1;
w = w_1;

%Initialization
X0  = X0_1;
X0m  = X0m_1;
w10 = w10_1;
w20 = w20_1;
Psi0 = Psi0_1;
xi0 = xi0_1;
Omega0 = Omega0_1;
rho0 = rho0_1;

init = [X0' X0m' w10' w20' Psi0' xi0' Omega0' rho0]';

%Adaptation gain and filter ctes
N = N_1;
c1 = c1_1;
c2 = c2_1;
d1 = d1_1;
d2 = d2_1;
Gamma = Gamma_2*eye(length(Psi0));
gamma = gamma_1;

%Simulation
[T_2,X_2] = ode23s('backstepping_red_direct',tfinal,init,'');
y_2       = X_2(:,1);
ym_2      = X_2(:,3);
Psi_2     = X_2(:,7:10);
rho_2     = X_2(:,end);
tilPsi_2  = -Psi_2 + Psis_2.*ones(length(Psi_2),length(Psis_2));
modPsi_2  = sqrt(sum(Psi_2.^2,2));
e0_2 =  y_2 - ym_2;

r_2 = 0;
for i=1:length(a)
    r_2 = r_2 + a(i)*sin(w(i).*T_2);
end

%Plot
run plot_bkst.m;

% %% -------------------------------------------------- Simulation 3 (planta)
% changed = 2;
% 
% %Plant
% kp = kp_2;
% a1 = a1_2;
% a0 = a0_2;
% P = tf(kp,[1 a1 a0]);
% 
% %Model
% km = km_1;
% am1 = am1_1;
% am0 = am0_1;
% Pm = tf(km,[1 am1 am0]);
% 
% %2DOF Control ideal parameters and lambda filter
% [t1, tn, t2, t2n, L] = find2DOFparameters(P,Pm,A0);
% Psis_2 = 1/t2n*[1 t1 tn t2];
% 
% %u and y filter
% ss_f = canon(ss(tf(1,L)), 'companion');
% Af = ss_f.A';
% Bf = ss_f.C';
% 
% % Reference
% dc = dc_1;
% a = a_1;
% w = w_1;
% 
% %Initialization
% X0  = X0_1;
% X0m  = X0m_1;
% w10 = w10_1;
% w20 = w20_1;
% Psi0 = Psi0_1;
% xi0 = xi0_1;
% Omega0 = Omega0_1;
% rho0 = rho0_1;
% 
% init = [X0' X0m' w10' w20' Psi0' xi0' Omega0' rho0]';
% 
% %Adaptation gain and filter ctes
% N = N_1;
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1*eye(length(Psi0));
% gamma = gamma_1;
% 
% %Simulation
% [T_2,X_2] = ode23s('backstepping_red_direct',tfinal,init,'');
% y_2       = X_2(:,1);
% ym_2      = X_2(:,3);
% Psi_2     = X_2(:,7:10);
% rho_2     = X_2(:,end);
% tilPsi_2  = -Psi_2 + Psis_2.*ones(length(Psi_2),length(Psis_2));
% modPsi_2  = sqrt(sum(Psi_2.^2,2));
% e0_2 =  y_2 - ym_2;
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% 
% %Plot
% run plot_bkst.m;
% 
% %% --------------------------------------------------- Simulation 4 (model)
% changed = 3;
% 
% %Plant
% kp = kp_1;
% a1 = a1_1;
% a0 = a0_1;
% P = tf(kp,[1 a1 a0]);
% 
% %Model
% km = km_2;
% am1 = am1_2;
% am0 = am0_2;
% Pm = tf(km,[1 am1 am0]);
% 
% %2DOF Control ideal parameters and lambda filter
% [t1, tn, t2, t2n, L] = find2DOFparameters(P,Pm,A0);
% Psis_2 = 1/t2n*[1 t1 tn t2];
% 
% %u and y filter
% ss_f = canon(ss(tf(1,L)), 'companion');
% Af = ss_f.A';
% Bf = ss_f.C';
% 
% % Reference
% dc = dc_1;
% a = a_1;
% w = w_1;
% 
% %Initialization
% X0  = X0_1;
% X0m  = X0m_1;
% w10 = w10_1;
% w20 = w20_1;
% Psi0 = Psi0_1;
% xi0 = xi0_1;
% Omega0 = Omega0_1;
% rho0 = rho0_1;
% 
% init = [X0' X0m' w10' w20' Psi0' xi0' Omega0' rho0]';
% 
% %Adaptation gain and filter ctes
% N = N_1;
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1*eye(length(Psi0));
% gamma = gamma_1;
% 
% %Simulation
% [T_2,X_2] = ode23s('backstepping_red_direct',tfinal,init,'');
% y_2       = X_2(:,1);
% ym_2      = X_2(:,3);
% Psi_2     = X_2(:,7:10);
% rho_2     = X_2(:,end);
% tilPsi_2  = -Psi_2 + Psis_2.*ones(length(Psi_2),length(Psis_2));
% modPsi_2  = sqrt(sum(Psi_2.^2,2));
% e0_2 =  y_2 - ym_2;
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% 
% %Plot
% run plot_bkst.m;
% 
% %% ------------------------------------------------------ Simulation 5 (y0)
% changed = 4;
% 
% %Plant
% kp = kp_1;
% a1 = a1_1;
% a0 = a0_1;
% P = tf(kp,[1 a1 a0]);
% 
% %Model
% km = km_1;
% am1 = am1_1;
% am0 = am0_1;
% Pm = tf(km,[1 am1 am0]);
% 
% %2DOF Control ideal parameters and lambda filter
% [t1, tn, t2, t2n, L] = find2DOFparameters(P,Pm,A0);
% Psis_2 = 1/t2n*[1 t1 tn t2];
% 
% %u and y filter
% ss_f = canon(ss(tf(1,L)), 'companion');
% Af = ss_f.A';
% Bf = ss_f.C';
% 
% % Reference
% dc = dc_1;
% a = a_1;
% w = w_1;
% 
% %Initialization
% X0  = X0_2;
% X0m  = X0m_1;
% w10 = w10_1;
% w20 = w20_1;
% Psi0 = Psi0_1;
% xi0 = xi0_1;
% Omega0 = Omega0_1;
% rho0 = rho0_1;
% 
% init = [X0' X0m' w10' w20' Psi0' xi0' Omega0' rho0]';
% 
% %Adaptation gain and filter ctes
% N = N_1;
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1*eye(length(Psi0));
% gamma = gamma_1;
% 
% %Simulation
% [T_2,X_2] = ode23s('backstepping_red_direct',tfinal,init,'');
% y_2       = X_2(:,1);
% ym_2      = X_2(:,3);
% Psi_2     = X_2(:,7:10);
% rho_2     = X_2(:,end);
% tilPsi_2  = -Psi_2 + Psis_2.*ones(length(Psi_2),length(Psis_2));
% modPsi_2  = sqrt(sum(Psi_2.^2,2));
% e0_2 =  y_2 - ym_2;
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% 
% %Plot
% run plot_bkst.m;
