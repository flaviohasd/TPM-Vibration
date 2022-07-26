function [theta1vg,vrelv,R1,delta,Tauv,Kte,Ktem,K12,C,Cm,TE,TEm,dTE] = energia(step,deltamax,F,double,alpha,pb,z1,z2,di1,di2,rb1,rb2,rp1,rp2,w,E,b,thetad,thetas,thetab1,thetab2,alpha10,alpha20,L,hf1,hf2,md1,md2,Ikbi,Ikai,Iksi,Kh,me,memod,Ra,v0)

    Kte = zeros(1,step);
    Ktem = zeros(1,step);
    C = zeros(1,step);
    Tauv = zeros(1,step);
    R1 = zeros(1,step);
    TEm = zeros(1,step);
    TE = zeros(1,step);
    dTE = zeros(1,step);
    Cm = zeros(1,step);
    delta = zeros(step,2);
    K12 = zeros(step,2);
    theta1vg = zeros(step,1);
    vrelv = zeros(step,1);

% Tabela para o calculo da rigidez da fundacao
    %         A              B            C             D          E       F
    X = [-5.574*1e-5,  -1.9986*1e-3, -2.3015*1e-4, 4.7702*1e-3,  0.0271 6.8045;   % L
         60.111*1e-5,  28.1*1e-3,    -83.431*1e-4, -9.9256*1e-3, 0.1624 0.9086;   % M
         -50.952*1e-5, 185.5*1e-3,   0.0538*1e-4,  53.3*1e-3,    0.2895 0.9236;   % P
         -6.2042*1e-5, 9.0889*1e-3,  -4.0964*1e-4, 7.8297*1e-3,  -0.1472 0.6904]; % Q

    XL1 = X(1,1)/thetab1^2 + X(1,2)*hf1^2 + X(1,3)*hf1/thetab1 + X(1,4)/thetab1 + X(1,5)*hf1 + X(1,6);
    XL2 = X(1,1)/thetab2^2 + X(1,2)*hf2^2 + X(1,3)*hf2/thetab2 + X(1,4)/thetab2 + X(1,5)*hf2 + X(1,6); 
    XM1 = X(2,1)/thetab1^2 + X(2,2)*hf1^2 + X(2,3)*hf1/thetab1 + X(2,4)/thetab1 + X(2,5)*hf1 + X(2,6);
    XM2 = X(2,1)/thetab2^2 + X(2,2)*hf2^2 + X(2,3)*hf2/thetab2 + X(2,4)/thetab2 + X(2,5)*hf2 + X(2,6); 
    XP1 = X(3,1)/thetab1^2 + X(3,2)*hf1^2 + X(3,3)*hf1/thetab1 + X(3,4)/thetab1 + X(3,5)*hf1 + X(3,6);
    XP2 = X(3,1)/thetab2^2 + X(3,2)*hf2^2 + X(3,3)*hf2/thetab2 + X(3,4)/thetab2 + X(3,5)*hf2 + X(3,6); 
    XQ1 = X(4,1)/thetab1^2 + X(4,2)*hf1^2 + X(4,3)*hf1/thetab1 + X(4,4)/thetab1 + X(4,5)*hf1 + X(4,6);
    XQ2 = X(4,1)/thetab2^2 + X(4,2)*hf2^2 + X(4,3)*hf2/thetab2 + X(4,4)/thetab2 + X(4,5)*hf2 + X(4,6);

    % LOOP PARA VARREDURA DO DESLOCAMENTO ANGULAR DO PINHAO
    thetat = @(t) w*t; % Deslocamento angular com relacao ao tempo
    i = 1;
    for t = linspace(0,(2*pi)/w,step) % Deslocamento angular do pinhao (OBS: usar step de no minimo 1e4 para evitar erros de iteracao da correcao do theta)
        theta1 = thetat(t); % Deslocamento angular no tempo especifico
        theta1vg(i) = rad2deg(theta1); % Vetor deslocamento angular (graus)
        
        % Ciclos dos dentes
            theta1corr = floor(theta1/(2*pi/z1))*(2*pi/z1); % theta1 corrigido para o ciclo do dente

        % Relacoes que definem os angulos e distancias de contato
            tdouble = [theta1corr; theta1corr+thetad]; % Vetor que define a fase de engrenamento duplo
            tsingle = [theta1corr+thetad; theta1corr+(thetad+thetas)]; % Vetor que define a fase de engrenamento simples

            theta = theta1-floor(theta1/(thetad+thetas))*(thetad+thetas); % theta1 corrigido para o ciclo do double + single
            thetadsd = theta1-floor(theta1/(2*thetad+thetas))*(2*thetad+thetas); % theta1 corrigido para o ciclo do double + single + double

            phi1 = theta+alpha10; % Angulo de engrenamento do (angulo que define o ponto de contato do par 1) (com relacao a engrenagem 1)
            phi1dsd = thetadsd+alpha10;  % phi para o calculo da vrel e coordenada generalizada
            phi1d = phi1 + 2*pi/z1; % phi1d = phi1+pb/rb1; % Angulo de engrenamento (angulo que define o ponto de contato do par 2) (com relacao a engrenagem 1)

            h1 = rb1*((thetab1+phi1)*cos(phi1)-sin(phi1)); % Largura do dente no ponto de contato (com relacao a engrenagem 1) (mm)
            h1d = rb1*((thetab1+phi1d)*cos(phi1d)-sin(phi1d));
            l1 = rb1*((thetab1+phi1)*sin(phi1)+cos(phi1)-cos(thetab1)); % Altura do dente no ponto de contato (com relacao a engrenagem 1) (mm)
            l1d = rb1*((thetab1+phi1d)*sin(phi1d)+cos(phi1d)-cos(thetab1));
            x1 = rb1*(thetab1+phi1); % Distancia do ponto de contato (do par 1) na linha de acao (com relacao a engrenagem 1) (mm)
            x1d = rb1*(thetab1+phi1d); % Distancia do ponto de contato (do par 2) na linha de acao (com relacao a engrenagem 1) (mm)

            x2 = L-x1; % Distancia do ponto de contato na linha de acao  (com relacao a engrenagem 2) (mm)
            x2d = L-x1d; % Distancia do ponto de contato (do par 2) na linha de acao (com relacao a engrenagem 2) (mm)
            phi2 = alpha20-(z1/z2)*theta; % phi2 = x2/rb2-thetab2; % Angulo de engrenamento (angulo que define o ponto de contato do par 1) (com relacao a engrenagem 2) (rad)
            phi2dsd = alpha20-(z1/z2)*thetadsd; % Angulo de engrenamento ajustado para o ciclo completo do dente (rad)
            phi2d = phi2-2*pi/z2; % phi2d = x2d/rb2-thetab2; % Angulo de engrenamento (angulo que define o ponto de contato do par 2) (com relacao a engrenagem 2) (rad)
            l2 = rb2*((thetab2+phi2)*sin(phi2)+cos(phi2)-cos(thetab2)); % (vulgo d) Altura do dente no ponto de contato (com relacao a engrenagem 1)  (mm)
            l2d = rb2*((thetab2+phi2d)*sin(phi2d)+cos(phi2d)-cos(thetab2));
            h2 = rb2*((thetab2+phi2)*cos(phi2)-sin(phi2)); % Largura do dente no ponto de contato (com relacao a engrenagem 2) (mm)
            h2d = rb2*((thetab2+phi2d)*cos(phi2d)-sin(phi2d));

            h1dsd = rb1*((thetab1+phi1dsd)*cos(phi1dsd)-sin(phi1dsd)); % h1 para calculo da coordenada generalizada (mm)
            l1dsd = rb1*((thetab1+phi1dsd)*sin(phi1dsd)+cos(phi1dsd)-cos(thetab1)); % l1 para calculo da coordenada generalizada (mm)
            x1dsd = rb1*(thetab1+phi1dsd); % x1 para o calculo da vrel (mm)
            x2dsd = L-x1dsd; % x1 para o calculo da vrel (mm)

        % Conversao da referencia do meio do dente para o ponto de contato
            alphax1 = phi1 + atan(h1/(rb1*cos(thetab1) + l1)); % Angulo do ponto de contato periodico do dente (rad)

            rx1 = rb1/cos(alphax1); % Distancia radial do primeiro ponto de contato (engrenagem 1) (mm)
            dx = rx1*sin(alphax1)-md1; % Distancia percorrida do engrenamento na linha de acao (mm)
            alphax1d = acos(sqrt((md1+dx)^2+rb1^2)*cos(alphax1)/(sqrt((md1+pb+dx)^2+rb1^2))); % Angulo do segundo ponto de contato (rad)
            rx1d = (rx1*sin(alphax1)+pb)/sin(alphax1d); % Distancia radial do segundo ponto de contato (engrenagem 1) (mm)

            rx2 = sqrt(rb2^2+(L-sqrt(rx1^2-rb1^2))^2); % Distancia radial do primeiro ponto de contato (engrenagem 2) (mm)
            alphax2 = acos(rb2/rx2); % (rad)
            rx2d = sqrt(rb2^2+(L-sqrt(rx1d^2-rb1^2))^2); % Distancia radial do segundo ponto de contato (engrenagem 2) (mm)

        % Coordenada Generalizada
            Tau =  (x1dsd/rb1)/tan(alpha) - 1;  % tan(alphax1dsd)/tan(alpha)-1;
            Tauv(i) = Tau;

        % Calculos dos desvios do dente devido a modificacao do perfil
            dx1 = x1-md1; % Distancia percorrida do ponto de contato 1 na linha de acao (util) com relacao ao pinhao (mm)
            dx2 = x2d-md2; % Distancia percorrida do ponto de contato 2 na linha de acao (util) com relacao a coroa (mm)

        % Rigidez pela energia
            uf1 = l1-h1*tan(phi1); Sf1 = thetab1*di1;
            uf2 = l2-h2*tan(phi2); Sf2 = thetab2*di2;
            uf1d = l1d-h1d*tan(phi1d); uf2d = l2d-h2d*tan(phi2d);
            Kf1 = (cos(phi1)^2/(E*b/1000))*(XL1*(uf1/Sf1)^2 + XM1*(uf1/Sf1) + XP1*(1+XQ1*tan(phi1)^2));
            Kf1d = (cos(phi1d)^2/(E*b/1000))*(XL1*(uf1d/Sf1)^2 + XM1*(uf1d/Sf1) + XP1*(1+XQ1*tan(phi1d)^2));
            Kf2 = (cos(phi2)^2/(E*b/1000))*(XL2*(uf2/Sf2)^2 + XM2*(uf2/Sf2) + XP2*(1+XQ2*tan(phi2)^2));
            Kf2d = (cos(phi2d)^2/(E*b/1000))*(XL2*(uf2d/Sf2)^2 + XM2*(uf2d/Sf2) + XP2*(1+XQ2*tan(phi2d)^2));

            if theta1 >= tdouble(1) && theta1 < tdouble(2) % Esta em fase 1 (duplo)
                K1 = Ikbi(phi1,thetab1)+Ikai(phi1,thetab1)+Iksi(phi1,thetab1)+Kf1;
                K2 = Ikbi(phi2,thetab2)+Ikai(phi2,thetab2)+Iksi(phi2,thetab2)+Kf2;
                K1d = Ikbi(phi1d,thetab1)+Ikai(phi1d,thetab1)+Iksi(phi1d,thetab1)+Kf1d;
                K2d = Ikbi(phi2d,thetab2)+Ikai(phi2d,thetab2)+Iksi(phi2d,thetab2)+Kf2d;
                Kte1 = 1/(1/Kh+K1+K2); Kte2 = 1/(1/Kh+K1d+K2d);
                Kte(i) = Kte1 + Kte2;
                K12(i,1) = Kte1; K12(i,2) = Kte2;

                deltax1 = (deltamax/1000)*(double-dx1)/double; % Erro do perfil devido modificacao com relacao a coroa (par 1) (mm)
                deltax2 = deltax1*(double-dx2)/(double-dx1); % Erro do perfil devido modificacao com relacao ao pinhao (par 2) (mm)
                delta(i,1) = deltax1/1000; delta(i,2) = deltax2/1000; % Erro dos pares engrenados (mm)

                Ktem(i) = F*Kte(i)/(F+Kte1*deltax1/1000+Kte2*deltax2/1000); % Rigidez equivalente  (N/m)
                
                R1(i) = Kte1/Kte(i);
            end

            if theta1 >= tsingle(1) && theta1 <= tsingle(2) % Esta em fase 2 (simples)
                K1 = Ikbi(phi1,thetab1)+Ikai(phi1,thetab1)+Iksi(phi1,thetab1)+Kf1;
                K2 = Ikbi(phi2,thetab2)+Ikai(phi2,thetab2)+Iksi(phi2,thetab2)+Kf2;
                Kte1 = 1/(1/Kh+K1+K2); Kte2 = 0;
                Kte(i) = Kte1 + Kte2;
                K12(i,1) = Kte1; K12(i,2) = Kte2;

                delta(i,1) = 0; delta(i,2) = 0;

                Ktem(i) = Kte(i); % Rigidez equivalente  (N/m)

                R1(i) = Kte1/Kte(i);
                mem = i;
            end

            if theta1 > (thetad+thetas) && theta1 <= (2*thetad+thetas) % Esta em fase 3 (duplo)
                R1(i) = 1- R1(i-mem);
            end

         % Erros de engrenamento
            % Energia
                TE(i) = F/(Kte1+Kte2); % Erro de transmissão teorico (m)
                dTE1 = (Kte1*deltax1/1000)/Kte(i); % Erro de transmissao estatica equivalente (par 1) (m)
                dTE2 = (Kte2*deltax2/1000)/Kte(i); % Erro de transmissao estatica equivalente (par 2) (m)
                dTE(i) = dTE1 + dTE2; % Erro de transmissao estatica equivalente total (m)
                TEm(i) = TE(i) + dTE(i); % Erro de transmissao estatica (apenas erros causados pelo TPM) (m)

        % Amortecimento
            Fn = F*cos(alphax1); % Forca normal a superficie do dente no ponto de contato (N)
            Ft = F*sin(alphax1); % Forca tangencial a superficie do dente no ponto de contato (N)
            vrel = abs(w*x1dsd/1000-(rp1/rp2)*w*x2dsd/1000); vrelv(i) = vrel;
            mu = 0.015; % 0.12*(Ft*Ra/(v0*vrel*rx1/1000))^0.25; % Coeficiente de atrito. Valor medio de mu: 0.015
            Fa = Fn*mu; % Forca de atrito
            B0 = Fn/Kte(i); % Amplitude da vibracao (equivalente a deformacao estatica);
            Wa = 4*Fa*B0; % Trabalho realizado pelo atrito
            U = (Fn^2)/(2*Kte(i)); % Energia perdida pelo atrito (equivalente à energia perdida na deformacao do dente);
            zeta = Wa/(4*pi*U); % Fator de amortecimento

            C(i) = 2*zeta*sqrt(me*Kte(i)); % Coeficiente de amortecimento (kg/(m²s)) ou (N*s/m³)
            Cm(i) = 2*zeta*sqrt(memod*Ktem(i)); % Coeficiente de amortecimento (kg/(m²s)) ou (N*s/m³)
            
        i = i + 1; % Contador para organizar os vetores
    end
end


