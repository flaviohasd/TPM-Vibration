function [dp1,de1,db1,di1,rd1,ra1,rb1,ri1,rp1,rf1,hf1,Sa1,Sd1,hr1,h1,hi1,dp2,de2,db2,di2,rd2,ra2,rb2,ri2,rp2,rf2,hf2,Sa2,Sd2,hr2,h2,hi2] = engrenagem(z1,z2,m,alpha)
x_shift = 0;
alphacrit = asin(sqrt(2/z1));   % Teste de angulo de pressao critico
if alphacrit > alpha
    alphacritdeg = alphacrit*180/pi;
    fprintf('Angulo de pressao menor que o critico: %.2f°! - Solução: Aumentar o numero de dentes. \n',alphacritdeg)
    return
end

e = 0.167*m; % Folga no pe do dente

% Engrenagem 1
    dp1 = m*z1; % Diametro primitivo 1 (mm)
    de1 = m*(z1+2); % Diametro externo 1 (mm)
    db1 = dp1*cos(alpha); % Diametro de base 1 (mm)
    di1 = m*(z1-2.334); % Diametro interno 1 (mm)

    rd1 = di1/2; % Raio do dedendo 1 (mm)
    rde1 = di1/2+e; % Raio do dedendo efetivo 1 (mm)
    ra1 = de1/2; % Raio do adendo 1 (mm)
    rb1 = db1/2; % Raio de base 1 (mm)
    ri1 = di1/2; % Raio interno 1 (mm)
    rp1 = dp1/2; % Raio primitivo 1 (mm)

    rf1 = di1*0.2/2; % Raio do furo (mm)
    hf1 = di1/(2*rf1); % Razao do furo pelo diametro interno

    sc1 = m*z1*sin((pi/2)/z1); % Espessura do dente no diametro primitivo 1 (mm)

    Sa1 = de1*(sc1/dp1+tan(alpha)-alpha-tan(acos(dp1*cos(alpha)/de1))+acos(dp1*cos(alpha)/de1)); % 1.74; % (mm)
    Sd1 = 2*rb1*sin((pi+4*x_shift*tan(alpha))/(2*z1) + tan(alpha) - alpha); % 5.8 no CAD (mm)
    hr1 = sqrt(rb1^2-(Sd1/2)^2) - sqrt(rd1^2-(Sd1/2)^2); % (mm)
    h1 = sqrt(ra1^2-(Sa1/2)^2) - sqrt(rd1^2-(Sd1/2)^2); % (mm)
    hi1 = (h1*Sd1-hr1*Sa1)/(Sd1-Sa1); % (mm)

% Engrenagem 2
    dp2 = m*z2; % Diametro primitivo 2 (mm)
    de2 = m*(z2+2); % Diametro externo 2 (mm)
    db2 = dp2*cos(alpha); % Diametro de base 2 (mm)
    di2 = m*(z2-2.334); % Diametro interno 2 (mm)

    rd2 = di2/2; % Raio do dedendo 2 (mm)
    rde2 = di2/2+e; % Raio do dedendo efetivo 2 (mm)
    ra2 = de2/2; % Raio do adendo 2 (mm)
    rb2 = db2/2; % Raio de base 2 (mm)
    ri2 = di2/2; % Raio interno 2 (mm)
    rp2 = dp2/2; % Raio primitivo 2 (mm)

    rf2 = di2*0.2/2; % Raio do furo (mm)
    hf2 = di2/(2*rf2); % Razao do furo pelo diametro interno

    sc2 = m*z2*sin((pi/2)/z2); % Espessura do dente no diametro primitivo 2 (mm)

    Sa2 = de2*(sc2/dp2+tan(alpha)-alpha-tan(acos(dp2*cos(alpha)/de2))+acos(dp2*cos(alpha)/de2)); % 1.83 (mm)
    alphaf = acos(rb2/rde2); % (rad)
    Sd2 = 2*rb2*sin((pi+4*x_shift*tan(alpha))/(2*z2) + tan(alpha) - alpha); % 6.03 no CAD (mm)
    hr2 = sqrt(rde2^2-(Sd2/2)^2) - sqrt(rd2^2-(Sd2/2)^2); % (mm)
    h2 = sqrt(ra2^2-(Sa2/2)^2) - sqrt(rd2^2-(Sd2/2)^2); % (mm)
    hi2 = (h2*Sd2-hr2*Sa2)/(Sd2-Sa2); % (mm)