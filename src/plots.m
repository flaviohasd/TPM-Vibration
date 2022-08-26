function [mm] = plots(Ampx,Ampa,RMSx,RMSa,Ampxtpm,Ampatpm,RMSxtpm,RMSatpm,Kte,deltamax,w,theta1vg,C,vrelv,thetad,thetas,R1,Kteint,KexData,KeyData,Cint,CxData,CyData,t1,t2,estados,estadostpm,tempo,acel,aceltpm,tvmf,tvmftpm,Tauv,Ktem,TE,TEm,step,z1,delta,ft,fttpm)

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
        
     % Ciclo do Tau
        pd = step/z1; % Quantidade de pontos em um ciclo de dente
        a = floor(pd*(z1-2+0.05)); b = floor(pd*(z1-0.4));
        R1mod=R1(floor(pd*5/100):b-a+floor(pd*5/100));
        Taur = linspace(-0.5,0.5,length(R1mod));
        Tau=linspace(-0.5,0.5,length(estados(a:b,1)));

    % RELACOES DE CONTATO
        figure('Name','Relacoes do Contato')
        subplot(3,1,1);
        plot(theta1vg,vrelv,'LineWidth',1.5);
        title('Velocidade relativa no ponto de contato')
        xlabel('\theta_1 (°)')
        ylabel('v_{rel} (m/s)')
        xlim([0 round((2*thetad+thetas)*180/pi,1)])

        subplot(3,1,2);
        plot(theta1vg,R1,'LineWidth',1.5);
        title('Razão de Transmissão')
        xlabel('\theta_1 (°)')
        ylabel('R')
        xlim([0 round((2*thetad+thetas)*180/pi,1)])
        
        subplot(3,1,3);
        plot(theta1vg,delta,'linewidth',1.5);
        title('\Delta')
        xlabel('\theta_1 (°)')
        ylabel('\Delta (m)')
        legend('\Delta_1','\Delta_2')
        xlim([0 round((2*thetad+thetas)*180/pi,1)])



    % RESULTADO DA INTERPOLACAO
        figure('Name','Resultado da Interpolação');
        subplot(2,1,1);
        plot( Kteint, KexData, KeyData );
        xlabel('t (s)')
        ylabel('k(t) (N/m)')
        legend('k(t)','k(t) interpolada')
        xlim([0 (2*pi)/w])

        subplot(2,1,2);
        plot( Cint, CxData, CyData );
        xlabel('t (s)')
        ylabel('c(t) (Nm/s)')
        legend('c(t)','c(t) interpolada')        
        xlim([0 (2*pi)/w])

    % SIMULACAO DO MODELO DINAMICO
        figure('Name','Dinâmica do engrenamento')
        subplot(2,2,1);
        plot(t1,estados(:,1),'LineWidth',1); hold on; plot(t2,estadostpm(:,1),'LineWidth',1);
        title({'Erro de transmissão dinâmico';'(Deslocamento relativo do engrenamento - na linha de ação)';''})
        xlabel('t (s)')
        ylabel('DTE (m)')
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
        ylabel('TVMF (N)')
        xlim([0 tempo])
        legend([h1(1) h2(1)],'\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

    % ANALISE EM COORDENADA NORMALIZADA
        figure('Name','Comparações com coordenada generalizada')
        subplot(3,2,1);
        plot(Tau,Kte(a:b),'-.','LineWidth',1.5); hold on;
        plot(Tau,Ktem(a:b),'LineWidth',1.5)
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
        xlabel('\Gamma')
        ylabel('k(t) (N/m)')
        xlim([-0.5 0.5])

        subplot(3,2,2);
        plot(Tau,TE(a:b),'-.','LineWidth',1.5); hold on;
        plot(Tau,TEm(a:b),'LineWidth',1.5)
        legend('TE',['TE_{M} ','(\Delta_{máx} =',' ',num2str(deltamax),' ',' (\mum))'])
        xlabel('\Gamma')
        ylabel('TE_{M} (m)')
        xlim([-0.5 0.5])

        subplot(3,2,3);
        plot(Tau,estados(a:b,1),'-.','LineWidth',1.5); hold on;
        plot(Tau,estadostpm(a:b,1),'LineWidth',1.5)
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
        xlabel('\Gamma')
        ylabel('DTE (m)')
        xlim([-0.5 0.5])

        subplot(3,2,4);
        plot(Tau,acel(a:b),'-.','LineWidth',1.5); hold on;
        plot(Tau,aceltpm(a:b),'LineWidth',1.5)
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
        xlabel('\Gamma')
        ylabel('a (m/s²)')
        xlim([-0.5 0.5])

        subplot(3,2,5);
        yyaxis right
        plot(Taur,R1mod,'k','LineWidth',1.5); hold on
        ylabel('R','color','k')
        ylim([0.05 1.2])
        yyaxis left
        plot(Taur,tvmf(a:b)'.*R1mod,'--','LineWidth',1.5);
        plot(Taur,tvmftpm(a:b)'.*R1mod,'color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
        legend('TVMF',['TVMF','(\Delta_{máx} =',' ',num2str(deltamax),')'],'Razão de Compartilhamento de Carga','Location','south','Orientation','horizontal')
        xlabel('\Gamma')
        ylabel('TVMF (N)')
        xlim([-0.5 0.5])
        ylim([0 max(max([tvmf(a:b)'.*R1mod,tvmftpm(a:b)'.*R1mod]))])
        ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
        
        subplot(3,2,6);
        plot(Tau,tvmf(a:b),'-.','LineWidth',1.5); hold on;
        plot(Tau,tvmftpm(a:b),'LineWidth',1.5)
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
        xlabel('\Gamma')
        ylabel('TVMF (N)')
        xlim([-0.5 0.5])

    % ANALISE DA FREQUENCIA
        figure('Name','Análise da frequência')
        plot(abs(ft(1:length(ft)/2))); hold on
        plot(abs(fttpm(1:length(fttpm)/2)));
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
        xlabel('f (Hz)')
        ylabel('Magnitude')

        mm = 0;

end
