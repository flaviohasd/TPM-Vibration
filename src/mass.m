function [me,memod] = mass(z1,z2,p,deltamax,cr,pb,b,rb1,rb2,rf1,rf2,ri1,ri2,Sd1,Sd2,hr1,hr2,hi1,hi2,h1,h2,Sa1,Sa2)

    ycg1 = ((Sd1*hr1^2)/2+(Sd1*(hi1-hr1)/2)*((hi1-hr1)/3+hr1)-(Sa1*(hi1-h1)/2)*((hi1-h1)/3+h1))/(Sd1*hr1+(Sd1*(hi1-hr1)/2)-(Sa1*(hi1-h1)/2)); % Tooth centroid 1 (mm)
    ycg2 = ((Sd2*hr2^2)/2+(Sd2*(hi2-hr2)/2)*((hi2-hr2)/3+hr2)-(Sa2*(hi2-h2)/2)*((hi2-h2)/3+h2))/(Sd2*hr2+(Sd2*(hi2-hr2)/2)-(Sa2*(hi2-h2)/2)); % Tooth centroid 2 (mm)
    
    ycg1mod = ((Sd1*hr1^2)/2+(Sd1*(hi1-hr1)/2)*((hi1-hr1)/3+hr1)-(Sa1*(hi1-h1)/2)*((hi1-h1)/3+h1)-(deltamax/1000)*(h1-hr1)*(2*(h1-hr1)/3+hr1))/(Sd1*hr1+(Sd1*(hi1-hr1)/2)-(Sa1*(hi1-h1)/2)-(deltamax/1000)*(h1-hr1)); % Modified tooth centroid 1 (mm)
    ycg2mod = ((Sd2*hr2^2)/2+(Sd2*(hi2-hr2)/2)*((hi2-hr2)/3+hr2)-(Sa2*(hi2-h2)/2)*((hi2-h2)/3+h2)-(deltamax/1000)*(h2-hr2)*(2*(h2-hr2)/3+hr2))/(Sd2*hr2+(Sd2*(hi2-hr2)/2)-(Sa2*(hi2-h2)/2)-(deltamax/1000)*(h1-hr1)); % Modified tooth centroid 2 (mm)
    
    A1 = Sd1*hr1 + Sd1*(hi1-hr1)/2 - Sa1*(hi1-h1)/2; % Transversal section area of tooth 1 (mm^2)
    A2 = Sd2*hr2 + Sd2*(hi2-hr2)/2 - Sa2*(hi2-h2)/2; % Transversal section area of tooth 2 (mm^2)
    Ar1 = (deltamax/1000)*(hi1-hr1); % Approximated area of the removed material of the modified tooth 1 (mm^2)
    Ar2 = (deltamax/1000)*(hi2-hr2); % Approximated area of the removed material of the modified tooth 2 (mm^2)
    
    V1 = A1*b/(1e9); % Volume of tooth 1 (m^3)
    V2 = A2*b/(1e9); % Volume of tooth 2 (m^3)
    Vr1 = Ar1*b/(1e9); % Approximated volume of the removed material of the modified tooth 1 (m^3)
    Vr2 = Ar2*b/(1e9); % Approximated volume of the removed material of the modified tooth 2 (m^3)
    
    mD1 = V1*p; % Mass of the tooth 1 (kg)
    mD2 = V2*p; % Mass of the tooth 2 (kg)
    mDr1 = Vr1*p; % Approximated mass of the removed material of the modified tooth (kg)
    mDr2 = Vr2*p; % Approximated mass of the removed material of the modified tooth (kg)

    mDt1 = mD1*z1; % Mass of all teeth of gear 1 (kg)
    mDt2 = mD2*z2; % Mass of all teeth of gear 2 (kg)

    mDt1mod = (mD1-mDr1)*z1; % Mass of all modified teeth of gear 1 (kg)
    mDt2mod = (mD2-mDr2)*z2; % Mass of all modified teeth of gear (kg)

    m1c = p*pi*(b/1000)*(ri1/1000)^2; m1cr = p*pi*(b/1000)*(rf1/1000)^2; % Mass of the gear body and mass removed from gear 1 (kg)
    m2c = p*pi*(b/1000)*(ri2/1000)^2; m2cr = p*pi*(b/1000)*(rf2/1000)^2; % Mass of the gear body and mass removed from gear 2 (kg)

    mt1 = mDt1 + m1c - m1cr; % Total mass of gear 1 (kg)
    mt2 = mDt2 + m2c - m2cr; % Total mass of gear 2 (kg)

    mt1mod = mDt1mod + m1c - m1cr; % Total mass of modified gear 1 (kg)
    mt2mod = mDt2mod + m2c - m2cr; % Total mass of modified gear (kg)

    Id1 = mDt1*(ycg1/1000+ri1/1000)^2; % Moment of inertia of the teeth from gear 1 (considering all mass at a singular point at the center of gravity)
    Id2 = mDt2*(ycg2/1000+ri2/1000)^2; % Moment of inertia of the teeth from gear 2 (considering all mass at a singular point at the center of gravity)

    Id1mod = mDt1mod*(ycg1mod/1000+ri1/1000)^2; % Moment of inertia of the teeth from modified gear 1 (considering all mass at a singular point at the center of gravity)
    Id2mod = mDt2mod*(ycg2mod/1000+ri2/1000)^2; % Moment of inertia of the teeth from modified gear 2 (considering all mass at a singular point at the center of gravity)

    Ic1 = m1c*(ri1/1000)^2/2 - m1cr*(rf1/1000)^2/2; % Moment of inertia of the gear body 1 (kg.m^2)
    Ic2 = m2c*(ri2/1000)^2/2 - m2cr*(rf2/1000)^2/2; % Moment of inertia of the gear body 2 (kg.m^2)

    J1 = Id1 + Ic1; % Polar moment of inertia of the gear 1 (kg.m^2)
    J2 = Id2 + Ic2; % Polar moment of inertia of the gear 2 (kg.m^2)

    J1mod = Id1mod + Ic1; % Polar moment of inertia of the modified gear 1 (kg.m^2)
    J2mod = Id2mod + Ic2; % Polar moment of inertia of the modified gear 2 (kg.m^2)

    me = J1*J2/(J2*(rb1/1000)^2 + J1*(rb2/1000)^2); % Equivalent mass of the system (kg)

    memod = J1mod*J2mod/(J2mod*(rb1/1000)^2 + J1mod*(rb2/1000)^2); % Equivalent mass of the modified system (kg)

end
