% function Ampxtpm = simulacao(deltamax)
    syms rx1 rx2 alpha1 alpha2 x

    %% PARAMETROS UTILIZADOS NO CODIGO
        tempo = 0.03; % Tempo de analise (s) (necessita ser ao menos 1 periodo do pinhao)
        step = 1e4; % Quantidade de iteracoes a serem realizadas no loop para calculo da rigidez e demais parametros (ATENCAO: utilizar no minimo 1e4 para evitar acumulo de erros das iteracoes do calculo da rigidez)
        deltamax = 18.5971; % Coeficiente de mofificacao do adendo maximo (µm)    %18.5886 para Fc / 18.5912 para acel
        x0 = [0; 0]; % Condicoes iniciais: [posicao(0) velocidade(0)] do modelo dinamico em espaco de estados

    %% CONSTANTES
        % GEOMETRIA E PROPRIEDADES DAS ENGRENAGENS
            z1 = 27; % N Dentes 1
            z2 = 35; % N Dentes 2
            m = 3; % Modulo (=25.4/P para P = [1/in])
            alpha = deg2rad(20); % Angulo de pressao (rad)
            b = 25; % Largura (mm)

            E = 206.8e9; % Modulo de elasticidade (Pa)
            v = 0.3; % Coeficiente de Poisson
            p = 7850; % Densidade do material (kg/m^3)
            Ra = 0.63/1e6; % Rugosidade da superficie (µm)

        % CONDICOES DE OPERACAO
            Pot = 80e3; % Potencia de entrada (W) (pinhao)
            W = 2000; % Rotacao do pinhao (rpm)

            p0 = 870;  % Densidade do lubrificante (kg/m^3)
            t0 = 60; % Temperatura de operacao (°C)
            v0 = 226.166; % Viscosidade cinematica (cst) -  320 (40°C) / 38.5 (100°C)

        % CALCULO DAS DIMENSOES:
            [dp1,de1,db1,di1,rd1,ra1,rb1,ri1,rp1,rf1,hf1,Sa1,Sd1,hr1,h1,hi1,dp2,de2,db2,di2,rd2,ra2,rb2,ri2,rp2,rf2,hf2,Sa2,Sd2,hr2,h2,hi2] = engrenagem(z1,z2,m,alpha);

        % CALCULO DAS FORCAS E VELOCIDADES
            [w,w1,w2,T,F,freq,ciclo] = operacao(Pot,W,db1,dp1,dp2);

    % RELACOES DO ENGRENAMENTO
       [pb,L,cr,PTH,md1,md2,thetad,thetas,alpha10,alpha20,thetab1,thetab2,double,single] = engrenamento(alpha,z1,z2,dp1,db1,ra1,rb1,dp2,ra2,rb2);

        f = 1;   % Fator compartilhamento de carga par 1 (engrenamento simples)
        f1 = 0.5*f; % Fator compartilhamento de carga par 1 (engrenamento duplo)
        f2 = f-f1; % Fator compartilhamento de carga par 2 (engrenamento duplo)

    %% MASSA E MOMENTO DE INERCIA
        [me,memod] = massa(z1,z2,p,deltamax,cr,pb,b,rb1,rb2,rf1,rf2,ri1,ri2,Sd1,Sd2,hr1,hr2,hi1,hi2,h1,h2,Sa1,Sa2);

    %% RIGIDEZ E AMORTECIMENTO
        [Ikbi,Iksi,Ikai,Kh] = integ(E,b,v);

        [theta1vg,vrelv,R1,delta,Tauv,Kte,Ktem,K12,C,Cm,TE,TEm,dTE] = energia(step,deltamax,F,double,alpha,pb,z1,z2,di1,di2,rb1,rb2,rp1,rp2,w,E,b,thetad,thetas,thetab1,thetab2,alpha10,alpha20,L,hf1,hf2,md1,md2,Ikbi,Ikai,Iksi,Kh,me,memod);

        % INTERPOLACAO DE K E C PARA ENTRADAS TEMPORAIS
            [tint,Tauint,Kteint,K1int,K2int,Ktmint,Cint,Cmint,Fint,d1int,d2int,KexData,KeyData,CxData,CyData] = interp(tempo,step,w,Tauv,Kte,K12,Ktem,C,Cm,R1,F,delta);

    %% MODELO DINAMICO
        % SISTEMA EM ESPACO DE ESTADOS
            [t1,t2,estados,estadostpm,acel,aceltpm,tvmf,tvmftpm] = modelo(step,tempo,Kteint,Ktmint,K1int,K2int,Cint,Cmint,d1int,d2int,F,me,memod,x0);

        % ANALISE DO SINAL
            % Regime permanente
                % Deslocamento
                    [RMSx,RMSxtpm,Ampx,Ampxtpm,Espectrox,Espectroxtpm] = sinal(estados,estadostpm,freq);
                    Fcx = Ampxtpm/RMSxtpm;

                % Aceleracao
                    [RMSa,RMSatpm,Ampa,Ampatpm,Espectroa,Espectroatpm] = sinal(acel,aceltpm,freq);
                                  

    %% PLOTS
        [mm] = plots(Ampx,Ampa,RMSx,RMSa,Ampxtpm,Ampatpm,RMSxtpm,RMSatpm,Kte,deltamax,w,theta1vg,C,vrelv,thetad,thetas,R1,Kteint,KexData,KeyData,Cint,CxData,CyData,t1,t2,estados,estadostpm,tempo,acel,aceltpm,tvmf,tvmftpm,Tauv,Ktem,TE,TEm);
        
% end
