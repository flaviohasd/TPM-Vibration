function [me,memod] = massa(z1,z2,p,deltamax,cr,pb,b,rb1,rb2,rf1,rf2,ri1,ri2,Sd1,Sd2,hr1,hr2,hi1,hi2,h1,h2,Sa1,Sa2)

    ycg1 = ((Sd1*hr1^2)/2+(Sd1*(hi1-hr1)/2)*((hi1-hr1)/3+hr1)-(Sa1*(hi1-h1)/2)*((hi1-h1)/3+h1))/(Sd1*hr1+(Sd1*(hi1-hr1)/2)-(Sa1*(hi1-h1)/2)); % Centroide do dente 1 (mm)
    ycg2 = ((Sd2*hr2^2)/2+(Sd2*(hi2-hr2)/2)*((hi2-hr2)/3+hr2)-(Sa2*(hi2-h2)/2)*((hi2-h2)/3+h2))/(Sd2*hr2+(Sd2*(hi2-hr2)/2)-(Sa2*(hi2-h2)/2)); % Centroide do dente 2 (mm)
    
    ycg1mod = ((Sd1*hr1^2)/2+(Sd1*(hi1-hr1)/2)*((hi1-hr1)/3+hr1)-(Sa1*(hi1-h1)/2)*((hi1-h1)/3+h1)-(deltamax/1000)*(h1-hr1)*(2*(h1-hr1)/3+hr1))/(Sd1*hr1+(Sd1*(hi1-hr1)/2)-(Sa1*(hi1-h1)/2)-(deltamax/1000)*(h1-hr1)); % Centroide do dente 1 modificado (mm)
    ycg2mod = ((Sd2*hr2^2)/2+(Sd2*(hi2-hr2)/2)*((hi2-hr2)/3+hr2)-(Sa2*(hi2-h2)/2)*((hi2-h2)/3+h2)-(deltamax/1000)*(h2-hr2)*(2*(h2-hr2)/3+hr2))/(Sd2*hr2+(Sd2*(hi2-hr2)/2)-(Sa2*(hi2-h2)/2)-(deltamax/1000)*(h1-hr1)); % Centroide do dente 2 modificado (mm)
    
    A1 = Sd1*hr1 + Sd1*(hi1-hr1)/2 - Sa1*(hi1-h1)/2; % Area da seccao transversal do dente 1 (mm^2)
    A2 = Sd2*hr2 + Sd2*(hi2-hr2)/2 - Sa2*(hi2-h2)/2; % Area da seccao transversal do dente 2 (mm^2)
    Ar = (deltamax/1000)*(cr-1)*pb/2; % Area aproximada do material removido na modificacao do perfil dos dentes (mm^2)

    V1 = A1*b/(1e9); % Volume do dente 1 (m^3)
    V2 = A2*b/(1e9); % Volume do dente 2 (m^3)
    Vr = Ar*b/(1e9); % Volume aproximado do material removido na modificacao do perfil dos dentes (m^3)

    mD1 = V1*p; % Massa do dente 1 (kg)
    mD2 = V2*p; % Massa do dente 2 (kg)
    mDr = Vr*p; % Massa aproximada que foi removida na modificacao do perfil dos dentes (kg)

    mDt1 = mD1*z1; % Massa de todos os dentes da engrenagem 1 (kg)
    mDt2 = mD2*z2; % Massa de todos os dentes da engrenagem 2 (kg)

    mDt1mod = (mD1-mDr)*z1; % Massa de todos os dentes da engrenagem 1 modificada (kg)
    mDt2mod = (mD2-mDr)*z2; % Massa de todos os dentes da engrenagem 2 modificada (kg)

    m1c = p*pi*(b/1000)*(ri1/1000)^2; m1cr = p*pi*(b/1000)*(rf1/1000)^2; % Massa do corpo e massa removida da engrenagem 1 (kg)
    m2c = p*pi*(b/1000)*(ri2/1000)^2; m2cr = p*pi*(b/1000)*(rf2/1000)^2; % Massa do corpo e massa removida da engrenagem 2 (kg)

    mt1 = mDt1 + m1c; % Massa total da engrenagem 1 (kg)
    mt2 = mDt2 + m2c; % Massa total da engrenagem 2 (kg)

    mt1mod = mDt1mod + m1c; % Massa total da engrenagem 1 modificada (kg)
    mt2mod = mDt2mod + m2c; % Massa total da engrenagem 2 modificada (kg)

    Id1 = mDt1*(ycg1/1000+ri1/1000)^2; % Momento de inercia dos dentes da engrenagem 1 (considerando toda massa concentrada na distancia do centro de gravidade - propriedade do raio de giracao)
    Id2 = mDt2*(ycg2/1000+ri2/1000)^2; % Momento de inercia dos dentes da engrenagem 2 (considerando toda massa concentrada na distancia do centro de gravidade - propriedade do raio de giracao)

    Id1mod = mDt1mod*(ycg1mod/1000+ri1/1000)^2; % Momento de inercia dos dentes da engrenagem 1 modificada (considerando toda massa concentrada na distancia do centro de gravidade - propriedade do raio de giracao)
    Id2mod = mDt2mod*(ycg2mod/1000+ri2/1000)^2; % Momento de inercia dos dentes da engrenagem 2 modificada (considerando toda massa concentrada na distancia do centro de gravidade - propriedade do raio de giracao)

    Ic1 = m1c*(ri1/1000)^2/2 - m1cr*(rf1/1000)^2/2; % Momento de inercia do corpo da engrenagem 1 (kg/m^2)
    Ic2 = m2c*(ri2/1000)^2/2 - m2cr*(rf2/1000)^2/2; % Momento de inercia do corpo da engrenagem 2 (kg/m^2)

    J1 = Id1 + Ic1; % Momento de inercia polar da engrenagem 1 (kg/m^2)
    J2 = Id2 + Ic2; % Momento de inercia polar da engrenagem 2 (kg/m^2)

    J1mod = Id1mod + Ic1; % Momento de inercia polar da engrenagem 1 modificada (kg/m^2)
    J2mod = Id2mod + Ic2; % Momento de inercia polar da engrenagem 2 modificada (kg/m^2)

    me = J1*J2/(J2*(rb1/1000)^2 + J1*(rb2/1000)^2); % Massa equivalente do sistema (kg/m^4)

    memod = J1mod*J2mod/(J2mod*(rb1/1000)^2 + J1mod*(rb2/1000)^2); % Massa equivalente do sistema modificado (kg/m^4)
end
