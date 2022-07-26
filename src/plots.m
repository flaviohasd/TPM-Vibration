function [mm] = plots(Ampx,Ampa,RMSx,RMSa,Ampxtpm,Ampatpm,RMSxtpm,RMSatpm,Kte,deltamax,w,theta1vg,C,vrelv,thetad,thetas,R1,Kteint,KexData,KeyData,Cint,CxData,CyData,t1,t2,estados,estadostpm,tempo,acel,aceltpm,tvmf,tvmftpm,Tauv,Ktem,TE,TEm)
    % SINAL
        disp(' ')
        disp('                    DESLOCAMENTO:')
        disp('(Sistema original: Delta = 0 micro m)')
        disp(string(['Amp = ', num2str(Ampx, '%10.5e\n'),' m,   RMS = ',num2str(RMSx, '%10.5e\n'),' m']))
        disp(string(['Fator de crista = ', num2str(Ampx/RMSx),',   Fator K = ',num2str(Ampx*RMSx, '%10.5e\n'),' m²']))
        disp(' ')
        disp(string(['(Sistema modificado: Delta = ', num2str(deltamax),' micro m)']))
        disp(string(['Amp = ', num2str(Ampxtpm, '%10.5e\n'),' m,   RMS = ',num2str(RMSxtpm, '%10.5e\n'),' m']))
        disp(string(['Fator de crista = ', num2str(Ampxtpm/RMSxtpm),',   Fator K = ',num2str(Ampxtpm*RMSxtpm, '%10.5e\n'),' m²']))
        disp(' ')
        disp('------------------------------------------------------')
        disp(' ')
        disp('                     ACELERAÇÃO:')
        disp('(Sistema original: Delta = 0 micro m)')
        disp(string(['Amp = ', num2str(Ampa, '%10.5e\n'),' m,   RMS = ',num2str(RMSa, '%10.5e\n'),' m']))
        disp(string(['Fator de crista = ', num2str(Ampa/RMSa),',   Fator K = ',num2str(Ampa*RMSa, '%10.5e\n'),' m²']))
        disp(' ')
        disp(string(['(Sistema modificado: Delta = ', num2str(deltamax),' micro m)']))
        disp(string(['Amp = ', num2str(Ampatpm, '%10.5e\n'),' m,   RMS = ',num2str(RMSatpm, '%10.5e\n'),' m']))
        disp(string(['Fator de crista = ', num2str(Ampatpm/RMSatpm),',   Fator K = ',num2str(Ampatpm*RMSatpm, '%10.5e\n'),' m²']))

        mm = 0;
        
% RIGIDEZ
    subplot(2,1,1);
    plot(theta1vg,Kte,'LineWidth',1.5);
    xlabel('\theta (°)')
    ylabel('K (N/m)')
    xlim([0 360])
    title('Rigidez - Energia')

    subplot(2,1,2);
    plot(theta1vg,C,'LineWidth',1.5);
    title('Amortecimento')
    xlabel('\theta (°)')
    ylabel('C (Nm/s)')
    xlim([0 360])

% RELACOES DE CONTATO
    figure('Name','Relacoes do Contato')
    subplot(2,1,1);
    plot(theta1vg,vrelv,'LineWidth',1.5);
    title('Velocidade relativa no ponto de contato')
    xlabel('\theta (°)')
    ylabel('V_{rel} (m/s)')
    xlim([0 round((2*thetad+thetas)*180/pi,1)])

    subplot(2,1,2);
    plot(theta1vg,R1,'LineWidth',1.5);
    title('Razão de Transmissibilidade')
    xlabel('\theta (°)')
    ylabel('R')
    xlim([0 round((2*thetad+thetas)*180/pi,1)])


% RESULTADO DA INTERPOLACAO
    figure('Name','Resultado da Interpolação');
    subplot(2,1,1);
    plot( Kteint, KexData, KeyData );
    xlabel('t (s)')
    ylabel('K (N/m)')
    legend('Kt vs. t','Kt interpolada')
    xlim([0 (2*pi)/w])

    subplot(2,1,2);
    plot( Cint, CxData, CyData );
    xlabel('t (s)')
    ylabel('C (Nm/s)')
    legend('C vs. t','C interpolada')        
    xlim([0 (2*pi)/w])

% SIMULACAO DO MODELO DINAMICO
    figure('Name','Dinâmica do engrenamento')
    subplot(2,2,1);
    plot(t1,estados(:,1),'LineWidth',1); hold on; plot(t2,estadostpm(:,1),'LineWidth',1);
    title({'Erro de transmissão dinâmico';'(Deslocamento relativo do engrenamento - na linha de ação)';''})
    xlabel('t (s)')
    ylabel('x (m)')
    xlim([0 tempo])
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

    subplot(2,2,2);
    plot(t1,acel,'LineWidth',1); hold on; plot(t2,aceltpm,'LineWidth',1);
    title({'Taxa de variação da valocidade do erro de transmissão';'(Aceleração do engrenamento - na linha de ação)'})
    xlabel('t (s)')
    ylabel('a (m/s²)')
    xlim([0 tempo])
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

    subplot(2,2,3);
    plot(t1,estados(:,2),'LineWidth',1); hold on; plot(t2,estadostpm(:,2),'LineWidth',1);
    title({'Taxa de variação do erro de transmissão';'(Velocidade do engrenamento - na linha de ação)'})
    xlabel('t (s)')
    ylabel('v (m/s)')
    xlim([0 tempo])
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

    subplot(2,2,4);
    h1 = plot(t1,tvmf,'color',[0 0.4470 0.7410],'LineWidth',1); hold on; h2 = plot(t2,tvmftpm,'colo',[0.8500 0.3250 0.0980],'LineWidth',1);
    title({'Força de engrenamento';'(Força atuante na linha de ação)'})
    xlabel('t (s)')
    ylabel('F (N)')
    xlim([0 tempo])
    legend([h1(1) h2(1)],'\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

% ANALISE EM COORDENADA NORMALIZADA
    figure('Name','Comparações com coordenada generalizada')
    subplot(3,2,1);
    plot(Tauv(5535:6115),Kte(5560:6140),'-.','LineWidth',1.5); hold on;
    plot(Tauv(5535:6115),Ktem(5560:6140),'LineWidth',1.5)
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
    xlabel('\Gamma')
    ylabel('Kt (N/m)')
    xlim([-0.5 0.5])

    subplot(3,2,2);
    plot(Tauv(5535:6115),TE(5560:6140),'-.','LineWidth',1.5); hold on;
    plot(Tauv(5535:6115),TEm(5560:6140),'LineWidth',1.5)
    legend('TE',['TEm ','(\Delta_{máx} =',' ',num2str(deltamax),' ',' (\mum))'])
    xlabel('\Gamma')
    ylabel('Erro de transmissão (m)')
    xlim([-0.5 0.5])

    subplot(3,2,3);
    plot(Tauv(5535:6115),estados(5560:6140,1),'-.','LineWidth',1.5); hold on;
    plot(Tauv(5535:6115),estadostpm(5560:6140,1),'LineWidth',1.5)
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
    xlabel('\Gamma')
    ylabel('x (m)')
    xlim([-0.5 0.5])

    subplot(3,2,4);
    plot(Tauv(5535:6115),acel(5560:6140),'-.','LineWidth',1.5); hold on;
    plot(Tauv(5535:6115),aceltpm(5560:6140),'LineWidth',1.5)
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
    xlabel('\Gamma')
    ylabel('a (m/s²)')
    xlim([-0.5 0.5])

    subplot(3,2,5);
    yyaxis right
    plot(Tauv(1:615),R1(1:615),'k','LineWidth',1.5); hold on
    ylabel('R','color','k')
    yyaxis left
    plot(Tauv(1:615),tvmf(5001:5615)'.*R1(1:615),'--','LineWidth',1.5);
    plot(Tauv(1:615),tvmftpm(5001:5615)'.*R1(1:615),'color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    legend('TVMF','TVMF (mod)','Razão de Transmissibilidade','Location','south','Orientation','horizontal')
    xlabel('\Gamma')
    ylabel('F (N)')
    xlim([-0.5 0.5])
    ylim([0 max(max([tvmf(5001:5615),tvmftpm(5001:5615)]))])

    subplot(3,2,6);
    plot(Tauv(5535:6115),tvmf(5560:6140),'-.','LineWidth',1.5); hold on;
    plot(Tauv(5535:6115),tvmftpm(5560:6140),'LineWidth',1.5)
    legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
    xlabel('\Gamma')
    ylabel('F (N)')
    xlim([-0.5 0.5])
end
