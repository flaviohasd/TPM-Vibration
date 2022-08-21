function [pb,L,cr,PTH,md1,md2,thetad,thetas,alpha10,alpha20,thetab1,thetab2,double,single] = engrenamento(alpha,z1,z2,dp1,db1,ra1,rb1,dp2,ra2,rb2)
    pb = pi*db1/z1; %Passo da base (igual para engrenagem 1 e 2) (mm)
    L = sqrt((dp1/2)^2-(rb1)^2) + sqrt((dp2/2)^2-(rb2)^2); % Distancia da reta tangente entre as duas circunferencias de base  (mm)
    cr = (sqrt(ra1^2-rb1^2) + sqrt(ra2^2-rb2^2) - (dp1+dp2)*sin(alpha)/2)/pb; % Razao de contato
    PTH = cr*pb; % Comprimento da linha de contato (mm)

    % Fases de engrenamento (duplo ou simples)
    double = (cr-1)*pb; % Comprimento na linha de contato onde havera dois pares de dentes engrenados (mm)
    single = pb - double; % Comprimento na linha de contato onde havera um par de dentes engrenados (mm)
    md1 = sqrt(ra1^2-rb1^2)-PTH; % Distancia nao engrenada da linha de acao 1 (mm)
    md2 = sqrt(ra2^2-rb2^2)-PTH; % Distancia nao engrenada da linha de acao 2 (mm)

    fase1 = atan(md1/rb1); % Duplo (rad)
    fase2 = atan((md1+double)/rb1); % Simples (rad)
    fase3 = atan((md1+double+single)/rb1); % Duplo (rad)
    alphar = fase3-fase1; % Range angular do ciclo de engrenamento do dente (rad)

    thetad = tan(acos(z1*cos(alpha)/(z1+2)))-2*pi/z1-tan(acos(z1*cos(alpha)/sqrt((z2+2)^2+(z1+z2)^2-2*(z2+2)*(z1+z2)*cos(acos(z2*cos(alpha)/(z2+2))-alpha)))); % Duracao angular da fase de duplo engrenamento (rad)
    thetas = 2*pi/z1-thetad; % Duracao angular da fase de engrenamento simples (rad)

    thetab1 = pi/(2*z1) + tan(alpha) - alpha; % Metade da distância angular da raiz do dente 1
    thetab2 = pi/(2*z2) + tan(alpha) - alpha; % Metade da distância angular da raiz do dente 2

    alpha10 = -pi/(2*z1)-tan(alpha)+alpha+tan(acos(z1*cos(alpha)/sqrt((z2+2)^2+(z1+z2)^2-2*(z2+2)*(z1+z2)*cos(acos(z2*cos(alpha)/(z2+2))-alpha)))); % Posicao angular inicial do dente 1 - par 1 (theta = 0)
    alpha20 = tan(acos(z2*cos(alpha)/(z2+2)))-pi/(2*z2)-tan(alpha)+alpha; % Posicao angular inicial do segundo dente 1 - par 2 (theta = 0)
end
