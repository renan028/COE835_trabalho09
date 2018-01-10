%-------- Print eps plots -----

close all;

%Set matlab interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%Set figures positions and size
fig_xpos = 500;
fig_ypos = 250;
fig_width = 600;
fig_height = 250;
fig_pos = [fig_xpos fig_ypos fig_width fig_height];
pos_pct = .075;

switch changed
    case 1
        str1 = strcat('$\Gamma=',num2str(Gamma_1),'$');
        str2 = strcat('$\Gamma=',num2str(Gamma_2),'$');
        file_name = strcat('Gamma',num2str(Gamma_1),'Gamma',num2str(Gamma_2));
        
    case 2 
        str1 = '$P_1$';
        str2 = '$P_2$';
        file_name = 'P1P2';
        
    case 3
        str1 = '$P_{m_1}$';
        str2 = '$P_{m_2}$';
        file_name = 'Pm1Pm2';
        
    case 4
        y0_str_1 = num2str(X0_1(1));
        y0_str_2 = num2str(X0_2(1));
        str1 = strcat('$y(0)=',y0_str_1,'$');
        str2 = strcat('$y(0)=',y0_str_2,'$');
        file_name = 'y01y02';
end

path_tilPsi = strcat('../relatorio/figs/tilPsi/',sim_str,file_name,'.eps');
path_modPsi = strcat('../relatorio/figs/modPsi/',sim_str,file_name,'.eps');
path_e0 = strcat('../relatorio/figs/e0/',sim_str,file_name,'.eps');
path_y = strcat('../relatorio/figs/y/',sim_str,file_name,'.eps');
path_rho = strcat('../relatorio/figs/rho/',sim_str,file_name,'.eps');

%--------------- Fig1: til_Psi -------------
figure(1);clf;
set(gcf,'position',[fig_pos(1:2) fig_pos(3) 2*fig_pos(4)]);

h1 = subplot(211);
plot(T_1,tilPsi_1);grid on;
title(strcat('$\tilde{\Psi}$ com~ ', str1));

h2 = subplot(212);
plot(T_2,tilPsi_2);grid on;
title(strcat('$\tilde{\Psi}$ com~ ', str2));
% h2.YLim = h1.YLim;

subplot(211);
legend('$\tilde{\Psi}_1$','$\tilde{\Psi}_2$','$\tilde{\Psi}_3$','$\tilde{\Psi}_4$','Location','SouthEast')
subplot(212);
legend('$\tilde{\Psi}_1$','$\tilde{\Psi}_2$','$\tilde{\Psi}_3$','$\tilde{\Psi}_4$','Location','SouthEast')

%Reduce gap btw subplots
set(h2,'Position',[h2.Position(1), h2.Position(2) + pos_pct*(h1.Position(2) - h2.Position(2)), h2.Position(3), h2.Position(4)]);

if PRINT
    print(path_tilPsi,'-depsc2','-painters')
end

%--------------- Fig2: mod Psi -------------
figure(2);clf;
set(gcf,'position',fig_pos);

plot(T_1,modPsi_1);grid on;hold on;
plot(T_2,modPsi_2);

if (changed == 2) || (changed == 3)
    plot(T_1,norm(Psis_1)*ones(1,length(T_1)));
    plot(T_2,norm(Psis_2)*ones(1,length(T_2)));
    hold off;
    legend(str1,str2,'$||\Psi_1^*||$','$||\Psi_2^*||$','Location','SouthEast');
else
    plot(T_1,norm(Psis_1)*ones(1,length(T_1)));hold off;
    legend(str1,str2,'$||\Psi^*||$','Location','SouthEast');
end

title('$||\Psi||$');

if PRINT
    print(path_modPsi,'-depsc2','-painters')
end

%--------------- Fig3: e -------------
figure(3);clf;
set(gcf,'position',fig_pos);

plot(T_1,e0_1);grid;hold on;
plot(T_2,e0_2);hold off;

title('$e_0$');
legend(str1,str2,'Location','SouthEast');

if PRINT
    print(path_e0,'-depsc2','-painters')
end

%--------------- Fig3: y -------------
figure(4);clf;
set(gcf,'position',[fig_pos(1:2) fig_pos(3) 2*fig_pos(4)]);

h1 = subplot(211);
plot(T_1,y_1,T_1,r_1);grid on;
title(strcat('$y$ com~ ', str1));

h2 = subplot(212);
plot(T_2,y_2,T_2,r_2);grid on;
title(strcat('$y$ com~ ', str2));
% h2.YLim = h1.YLim;

subplot(211);
legend('$y$','$r$','Location','SouthEast')
subplot(212);
legend('$y$','$r$','Location','SouthEast')

set(gcf,'position',fig_pos);

if PRINT
    print(path_y,'-depsc2','-painters')
end

%--------------- Fig4: rho -------------
figure(3);clf;
set(gcf,'position',fig_pos);

plot(T_1,rho_1);grid;hold on;
plot(T_2,rho_2);

if (changed == 2) || (changed == 3)
    plot(T_1,1/Psis_1(1)*ones(1,length(T_1)));
    plot(T_2,1/Psis_2(1)*ones(1,length(T_2)));
    hold off;
    legend(str1,str2,'$\rho_1^*$','$\rho_2^*$','Location','SouthEast');
else
    plot(T_1,1/Psis_1(1)*ones(1,length(T_1)));hold off;
    legend(str1,str2,'$\rho^*$','Location','SouthEast');
end

title('$\rho$');

if PRINT
    print(path_rho,'-depsc2','-painters')
end
