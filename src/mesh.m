function [pb,L,cr,PTH,md1,md2,thetad,thetas,alpha10,alpha20,thetab1,thetab2,double,single] = mesh(alpha,z1,z2,dp1,db1,ra1,rb1,dp2,ra2,rb2)
    pb = pi*db1/z1; % Base step (equal to gear 1 and 2) (mm)
    L = sqrt((dp1/2)^2-(rb1)^2) + sqrt((dp2/2)^2-(rb2)^2); % Distance of the tangent line between the two base circunferences  (mm)
    cr = (sqrt(ra1^2-rb1^2) + sqrt(ra2^2-rb2^2) - (dp1+dp2)*sin(alpha)/2)/pb; % Contact ratio
    PTH = cr*pb; % Length of the contact line (mm)

    % Meshing contact phases (duplo ou simples)
    double = (cr-1)*pb; % Length of the contact line of the double contact phase (mm)
    single = pb - double; % Length of the line of contact of the single contact phase (mm)
    md1 = sqrt(ra1^2-rb1^2)-PTH; % Distance not in mesh in the line of action 1 (mm)
    md2 = sqrt(ra2^2-rb2^2)-PTH; % Distance not in mesh in the line of action 2 (mm)

    fase1 = atan(md1/rb1); % Double (rad)
    fase2 = atan((md1+double)/rb1); % Single (rad)
    fase3 = atan((md1+double+single)/rb1); % Double (rad)
    alphar = fase3-fase1; % Angular range of the tooth contact cycle (rad)

    thetad = tan(acos(z1*cos(alpha)/(z1+2)))-2*pi/z1-tan(acos(z1*cos(alpha)/sqrt((z2+2)^2+(z1+z2)^2-2*(z2+2)*(z1+z2)*cos(acos(z2*cos(alpha)/(z2+2))-alpha)))); % Angular length of the double contact phase (rad)
    thetas = 2*pi/z1-thetad; % Angular length of the single contact phase (rad)

    thetab1 = pi/(2*z1) + tan(alpha) - alpha; % Half angular distance of the tooth root 1
    thetab2 = pi/(2*z2) + tan(alpha) - alpha; % Half angular distance of the tooth root 2

    alpha10 = -pi/(2*z1)-tan(alpha)+alpha+tan(acos(z1*cos(alpha)/sqrt((z2+2)^2+(z1+z2)^2-2*(z2+2)*(z1+z2)*cos(acos(z2*cos(alpha)/(z2+2))-alpha)))); % Initial angular position of the tooth 1 - pair 1 (theta = 0)
    alpha20 = tan(acos(z2*cos(alpha)/(z2+2)))-pi/(2*z2)-tan(alpha)+alpha; % Initial angular position of the tooth 1 - pair 2 (theta = 0)
end
