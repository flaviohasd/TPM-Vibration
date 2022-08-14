function [Ikbi,Iksi,Ikai,Kh] = integ(E,b,v)
    syms alpha1 alpha2 x

    Fkbi = @(x) 3*((1+cos(alpha1)*((alpha2-x)*sin(x)-cos(x)))^2*(alpha2-x)*cos(x))/(2*E*(b/1000)*(sin(x)+(alpha2-x)*cos(x))^3); % 1/kb = integral(FUNCAO)
    iFkbi = @(x)int(Fkbi(x),-alpha1,alpha2); % 1/kb = INTEGRAL(funcao)
    Ikbi = matlabFunction(iFkbi(x)); % 1/KB = integral(funcao) (Rigidez de flexao)

    Fksi = @(x) (1.2*(1+v)*(alpha2-x)*cos(x)*cos(alpha1)^2)/(E*(b/1000)*(sin(x)+(alpha2-x)*cos(x))); % 1/ks = integral(FUNCAO)
    iFksi = @(x)int(Fksi(x),-alpha1,alpha2); % 1/ks = INTEGRAL(funcao)
    Iksi = matlabFunction(iFksi(x));  % 1/KS = integral(funcao) (Rigidez de compressao axial)

    Fkai = @(x) ((alpha2-x)*cos(x)*sin(alpha1)^2)/(2*E*(b/1000)*(sin(x)+(alpha2-x)*cos(x))); % 1/ka = integral(FUNCAO)
    iFkai = @(x)int(Fkai(x),-alpha1,alpha2); % 1/ka = INTEGRAL(funcao)
    Ikai = matlabFunction(iFkai(x));  % 1/KA = integral(funcao) (Rigidez de cisalhamento)

    Kh = pi*E*(b/1000)/(4*(1-v^2));  % Rigidez de contato Hertziano
end
