function [Ikbi,Iksi,Ikai,Kh] = integr(E,b,v)
    syms alpha1 alpha2 x

    Fkbi = @(x) 3*((1+cos(alpha1)*((alpha2-x)*sin(x)-cos(x)))^2*(alpha2-x)*cos(x))/(2*E*(b/1000)*(sin(x)+(alpha2-x)*cos(x))^3); % 1/kb = integral(FUNCTION)
    iFkbi = @(x)int(Fkbi(x),-alpha1,alpha2); % 1/kb = INTEGRAL(function)
    Ikbi = matlabFunction(iFkbi(x)); % 1/KB = integral(function) (Bending stiffness)

    Fksi = @(x) (1.2*(1+v)*(alpha2-x)*cos(x)*cos(alpha1)^2)/(E*(b/1000)*(sin(x)+(alpha2-x)*cos(x))); % 1/ks = integral(FUNCTION)
    iFksi = @(x)int(Fksi(x),-alpha1,alpha2); % 1/ks = INTEGRAL(function)
    Iksi = matlabFunction(iFksi(x));  % 1/KS = integral(function) (Shearing stiffness)

    Fkai = @(x) ((alpha2-x)*cos(x)*sin(alpha1)^2)/(2*E*(b/1000)*(sin(x)+(alpha2-x)*cos(x))); % 1/ka = integral(FUNCTION)
    iFkai = @(x)int(Fkai(x),-alpha1,alpha2); % 1/ka = INTEGRAL(function)
    Ikai = matlabFunction(iFkai(x));  % 1/KA = integral(function) (Axial compression stiffness)

    Kh = pi*E*(b/1000)/(4*(1-v^2));  % Hertzian contact stiffness
end
