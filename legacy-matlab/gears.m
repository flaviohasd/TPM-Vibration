function [dp1,de1,db1,di1,rd1,ra1,rb1,ri1,rp1,rf1,hf1,Sa1,Sd1,hr1,h1,hi1,dp2,de2,db2,di2,rd2,ra2,rb2,ri2,rp2,rf2,hf2,Sa2,Sd2,hr2,h2,hi2] = gears(z1,z2,m,alpha)
x_shift = 0;
alphacrit = asin(sqrt(2/z1));   % Teste de angulo de pressao critico
if alphacrit > alpha
    alphacritdeg = alphacrit*180/pi;
    fprintf('Pressure angle lower than critic: %.2f°! - Solution: Increase the number of teeth. \n',alphacritdeg)
    return
end

e = 0.167*m; % Clearance in the tooth root

% Gear 1
    dp1 = m*z1; % Primitive diameter 1 (mm)
    de1 = m*(z1+2); % External diameter 1 (mm)
    db1 = dp1*cos(alpha); % Base diameter 1 (mm)
    di1 = m*(z1-2.334); % Internal diameter 1 (mm)

    rd1 = di1/2; % Dedendum radius 1 (mm)
    rde1 = di1/2+e; % Effective dedendum radius 1 (mm)
    ra1 = de1/2; % Addendum radius 1 (mm)
    rb1 = db1/2; % Base radius 1 (mm)
    ri1 = di1/2; % Internal radius 1 (mm)
    rp1 = dp1/2; % Primitive radius 1 (mm)

    rf1 = di1*0.2/2; % Bore radius (mm) (set as 20% of the internal radius)
    hf1 = di1/(2*rf1); % Ratio of the internal diameter and bore diameter

    sc1 = m*z1*sin((pi/2)/z1); % Tooth width at the primitive point 1 (mm)

    Sa1 = de1*(sc1/dp1+tan(alpha)-alpha-tan(acos(dp1*cos(alpha)/de1))+acos(dp1*cos(alpha)/de1)); % 1.74; % (mm)
    Sd1 = 2*rb1*sin((pi+4*x_shift*tan(alpha))/(2*z1) + tan(alpha) - alpha); % (mm)
    hr1 = sqrt(rb1^2-(Sd1/2)^2) - sqrt(rd1^2-(Sd1/2)^2); % (mm)
    h1 = sqrt(ra1^2-(Sa1/2)^2) - sqrt(rd1^2-(Sd1/2)^2); % (mm)
    hi1 = (h1*Sd1-hr1*Sa1)/(Sd1-Sa1); % (mm)

% Gear 2
    dp2 = m*z2; % Primitive diameter 2 (mm)
    de2 = m*(z2+2); % External diameter 2 (mm)
    db2 = dp2*cos(alpha); % Base diameter 2 (mm)
    di2 = m*(z2-2.334); % Internal diameter 2 (mm)

    rd2 = di2/2; % Dedendum radius 2 (mm)
    rde2 = di2/2+e; % Effective dedendum radius 2 (mm)
    ra2 = de2/2; % Addendum radius 2 (mm)
    rb2 = db2/2; % Base radius 2 (mm)
    ri2 = di2/2; % Internal radius 2 (mm)
    rp2 = dp2/2; % Primitive radius 2 (mm)

    rf2 = di2*0.2/2; % Bore radius (mm) (set as 20% of the internal radius)
    hf2 = di2/(2*rf2); % Ratio of the internal diameter and bore diameter

    sc2 = m*z2*sin((pi/2)/z2); % Tooth width at the primitive point 2 (mm)

    Sa2 = de2*(sc2/dp2+tan(alpha)-alpha-tan(acos(dp2*cos(alpha)/de2))+acos(dp2*cos(alpha)/de2)); % (mm)
    alphaf = acos(rb2/rde2); % (rad)
    Sd2 = 2*rb2*sin((pi+4*x_shift*tan(alpha))/(2*z2) + tan(alpha) - alpha); % CAD (mm)
    hr2 = sqrt(rde2^2-(Sd2/2)^2) - sqrt(rd2^2-(Sd2/2)^2); % (mm)
    h2 = sqrt(ra2^2-(Sa2/2)^2) - sqrt(rd2^2-(Sd2/2)^2); % (mm)
    hi2 = (h2*Sd2-hr2*Sa2)/(Sd2-Sa2); % (mm)