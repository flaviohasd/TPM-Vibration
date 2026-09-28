function [mm] = plots(Ampx,Ampa,RMSx,RMSa,Ampxtpm,Ampatpm,RMSxtpm,RMSatpm,Kte,deltamax,w,theta1vg,C,vrelv,thetad,thetas,R1,Kteint,KexData,KeyData,Cint,CxData,CyData,t,estados,estadostpm,tempo,acel,aceltpm,tvmf,tvmftpm,Tauv,Ktem,TE,TEm,step,z1,delta,ft,fttpm)

    % SINAL
        disp(' ')
        disp('                    DISPLACEMENT:')
        disp('(Original system: Delta = 0 micro m)')
        disp(string(['Amp = ', num2str(Ampx, '%10.5e\n'),' m,   RMS = ',num2str(RMSx, '%10.5e\n'),' m']))
        disp(string(['Crest factor = ', num2str(Ampx/RMSx),',   K factor = ',num2str(Ampx*RMSx, '%10.5e\n'),' m²']))
        disp(' ')
        disp(string(['(Modified system: Delta = ', num2str(deltamax),' micro m)']))
        disp(string(['Amp = ', num2str(Ampxtpm, '%10.5e\n'),' m,   RMS = ',num2str(RMSxtpm, '%10.5e\n'),' m']))
        disp(string(['Crest factor = ', num2str(Ampxtpm/RMSxtpm),',   K factor = ',num2str(Ampxtpm*RMSxtpm, '%10.5e\n'),' m²']))
        disp(' ')
        disp('------------------------------------------------------')
        disp(' ')
        disp('                     ACCELERATION:')
        disp('(Original system: Delta = 0 micro m)')
        disp(string(['Amp = ', num2str(Ampa, '%10.5e\n'),' m,   RMS = ',num2str(RMSa, '%10.5e\n'),' m']))
        disp(string(['Crest factor = ', num2str(Ampa/RMSa),',   K factor = ',num2str(Ampa*RMSa, '%10.5e\n'),' m²']))
        disp(' ')
        disp(string(['(Modified System: Delta = ', num2str(deltamax),' micro m)']))
        disp(string(['Amp = ', num2str(Ampatpm, '%10.5e\n'),' m,   RMS = ',num2str(RMSatpm, '%10.5e\n'),' m']))
        disp(string(['Crest factor = ', num2str(Ampatpm/RMSatpm),',   K factor = ',num2str(Ampatpm*RMSatpm, '%10.5e\n'),' m²']))
        
     % NORMALISED COORDINATE
        pd = step/z1; % Quantity of points in a tooth cycle
        a = floor(pd*(z1-2+0.05)); b = floor(pd*(z1-0.4));
        R1mod=R1(floor(pd*5/100):b-a+floor(pd*5/100));
        Taur = linspace(-0.5,0.5,length(R1mod));
        Tau=linspace(-0.5,0.5,length(estados(a:b,1)));

    % CONTACT RELATIONS
        figure('Name','Contact relations')
        subplot(3,1,1);
        plot(theta1vg,vrelv,'LineWidth',1.5);
        title('Relative velocity at the contact point')
        xlabel('\theta_1 (°)')
        ylabel('v_{rel} (m/s)')
        xlim([0 round((2*thetad+thetas)*180/pi,1)])

        subplot(3,1,2);
        plot(theta1vg,R1,'LineWidth',1.5);
        title('Load sharing ratio')
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



    % INTERPOLATION RESULTS
        figure('Name','Interpolation results');
        subplot(2,1,1);
        plot( Kteint, KexData, KeyData );
        xlabel('t (s)')
        ylabel('k(t) (N/m)')
        legend('k(t)','k(t) interpolated')
        xlim([0 (2*pi)/w])

        subplot(2,1,2);
        plot( Cint, CxData, CyData );
        xlabel('t (s)')
        ylabel('c(t) (Nm/s)')
        legend('c(t)','c(t) interpolated')        
        xlim([0 (2*pi)/w])

    % SIMULATION OF THE DYNAMIC MODEL
        figure('Name','Meshing dynamics')
        subplot(2,2,1);
        plot(t,estados(:,1),'LineWidth',1); hold on; plot(t,estadostpm(:,1),'LineWidth',1);
        title({'Dynamic transmission error';'(Relative meshing displacement at the line of action)';''})
        xlabel('t (s)')
        ylabel('DTE (m)')
        xlim([0 tempo])
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
                
        subplot(2,2,2);
        plot(t,acel,'LineWidth',1); hold on; plot(t,aceltpm,'LineWidth',1);
        title({'Rate of change of the dynamic transmission error velocity';'(Meshing acceleration at the line of action)'})
        xlabel('t (s)')
        ylabel('a (m/s²)')
        xlim([0 tempo])
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

        subplot(2,2,3);
        plot(t,estados(:,2),'LineWidth',1); hold on; plot(t,estadostpm(:,2),'LineWidth',1);
        title({'Rate of change of the dynamic transmission error';'(Meshing velocity at the line of action)'})
        xlabel('t (s)')
        ylabel('v (m/s)')
        xlim([0 tempo])
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

        subplot(2,2,4);
        h1 = plot(t,tvmf,'color',[0 0.4470 0.7410],'LineWidth',1); hold on; h2 = plot(t,tvmftpm,'colo',[0.8500 0.3250 0.0980],'LineWidth',1);
        title({'Dynamic meshing force';'(Force acting in the line of action)'})
        xlabel('t (s)')
        ylabel('DMF (N)')
        xlim([0 tempo])
        legend([h1(1) h2(1)],'\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')

    % ANALYSIS WITH NORMALISED COORDINATE
        figure('Name','Comparisons with normalised coordinate')
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

    % FREQUENCY ANALYSIS
        figure('Name','Frequency analysis')
        plot(abs(ft(1:length(ft)/2))); hold on
        plot(abs(fttpm(1:length(fttpm)/2)));
        legend('\Delta_{máx} = 0 (\mum)',['\Delta_{máx} =',' ',num2str(deltamax),' ','(\mum)'],'Location','north','Orientation','horizontal')
        xlabel('f (Hz)')
        ylabel('Magnitude')

        mm = 0;

end
