function [theta1vg,vrelv,R1,delta,Tauv,Kte,Ktem,K12,C,Cm,TE,TEm,dTE] = energy(step,deltamax,F,double,alpha,pb,z1,z2,di1,di2,rb1,rb2,rp1,rp2,w,E,b,thetad,thetas,thetab1,thetab2,alpha10,alpha20,L,hf1,hf2,md1,md2,Ikbi,Ikai,Iksi,Kh,me,memod,Ra,v0)

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

% TABLE FOR FOUNDATION STIFFNESS CALCULATION
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

    % LOOP THROUGH THE ANGULAR DISPLACEMENT OF THE PINION
    thetat = @(t) w*t; % Angular displacement in time
    i = 1;
    for t = linspace(0,(2*pi)/w,step) % Angular displacement of the pinion
        theta1 = thetat(t); % Angular displacement at a specific time
        theta1vg(i) = rad2deg(theta1); % Vector of angular displacement (degree)
        
        % Tooth cycle
            theta1corr = floor(theta1/(2*pi/z1))*(2*pi/z1); % theta1 corrected for the tooth cycle

        % Relations that define the angles and the contact distances
            tdouble = [theta1corr; theta1corr+thetad]; % Vector that defines the double contact phase
            tsingle = [theta1corr+thetad; theta1corr+(thetad+thetas)]; % Vector that defines the single contact phase

            theta = theta1-floor(theta1/(thetad+thetas))*(thetad+thetas); % theta1 corrected for the double + single cycle
            thetadsd = theta1-floor(theta1/(2*thetad+thetas))*(2*thetad+thetas); % theta1 corrected for the double + single + double cycle

            phi1 = theta+alpha10; % Meshing angle (defines the contact point of pair of teeth 1) (with relation to gear 1)
            phi1dsd = thetadsd+alpha10;  % phi for vrel calculation and normalised coordinate
            phi1d = phi1 + 2*pi/z1; % phi1d = phi1+pb/rb1; % Meshing angle (defines the contact point of pair of teeth 2) (with relation to gear 1)

            h1 = rb1*((thetab1+phi1)*cos(phi1)-sin(phi1)); % Tooth width at the contact point (with relation to gear 1) (mm)
            h1d = rb1*((thetab1+phi1d)*cos(phi1d)-sin(phi1d));
            l1 = rb1*((thetab1+phi1)*sin(phi1)+cos(phi1)-cos(thetab1)); % Tooth hight at the contact point (with relation to gear 1) (mm)
            l1d = rb1*((thetab1+phi1d)*sin(phi1d)+cos(phi1d)-cos(thetab1));
            x1 = rb1*(thetab1+phi1); % Distance of the contact point (of pair 1) in the line of action (with relation to gear 1) (mm)
            x1d = rb1*(thetab1+phi1d); % Distance of the contact point (of pair 2) in the line of action (with relation to gear 1) (mm)

            x2 = L-x1; % Distance of the contact point in the line of action (with relation to gear 2) (mm)
            x2d = L-x1d; % Distance of the contact point in the line of action (of pair 2) in the line of action (with relation to gear 2) (mm)
            phi2 = alpha20-(z1/z2)*theta; % phi2 = x2/rb2-thetab2; % Meshing angle (defines the contact point of the pair of teeth 1) (with relation to gear2) (rad)
            phi2dsd = alpha20-(z1/z2)*thetadsd; % Meshing angle ajusted to the complete cycle of contact of the teeth (rad)
            phi2d = phi2-2*pi/z2; % phi2d = x2d/rb2-thetab2; % Meshing angle (defines the contact point of the pair of teeth 2) (with relation to gear 2) (rad)
            l2 = rb2*((thetab2+phi2)*sin(phi2)+cos(phi2)-cos(thetab2)); % (aka d) Tooth hight at the contact point (with relation to gear 1)  (mm)
            l2d = rb2*((thetab2+phi2d)*sin(phi2d)+cos(phi2d)-cos(thetab2));
            h2 = rb2*((thetab2+phi2)*cos(phi2)-sin(phi2)); % Tooth width at the contact point (with relation to gear 2) (mm)
            h2d = rb2*((thetab2+phi2d)*cos(phi2d)-sin(phi2d));

            h1dsd = rb1*((thetab1+phi1dsd)*cos(phi1dsd)-sin(phi1dsd)); % h1 for the calculation of the normalised coordinate (mm)
            l1dsd = rb1*((thetab1+phi1dsd)*sin(phi1dsd)+cos(phi1dsd)-cos(thetab1)); % l1 for the calculation of the normalised coordinate (mm)
            x1dsd = rb1*(thetab1+phi1dsd); % x1 for vrel calculation (mm)
            x2dsd = L-x1dsd; % x2 for vrel calculation (mm)

        % Converting the reference from the middle of the tooth to the contact point
            alphax1 = phi1 + atan(h1/(rb1*cos(thetab1) + l1)); % Angle of the periodic tooth contact point (rad)
            alphax1dsd = phi1dsd + atan(h1dsd/(rb1*cos(thetab1) + l1dsd));

            rx1 = rb1/cos(alphax1); % Radial distance of the first point of contact (gear 1) (mm)
            rx1dsd = rb1/cos(alphax1dsd);
            dx = rx1*sin(alphax1)-md1; % Distance traveled from the gear to the line of action (mm)
            alphax1d = acos(sqrt((md1+dx)^2+rb1^2)*cos(alphax1)/(sqrt((md1+pb+dx)^2+rb1^2))); % Angle of the second contact point (rad)
            rx1d = (rx1*sin(alphax1)+pb)/sin(alphax1d); % Radial distance of the second contact point (gear 1) (mm)

            rx2 = sqrt(rb2^2+(L-sqrt(rx1^2-rb1^2))^2); % Radial distance from the first contact point (gear 2) (mm)
            alphax2 = acos(rb2/rx2); % (rad)
            rx2d = sqrt(rb2^2+(L-sqrt(rx1d^2-rb1^2))^2); % Radial distance of the second contact point (gear 2) (mm)

        % Normalized Coordinate
            Tau =  (x1dsd/rb1)/tan(alpha) - 1;  % tan(alphax1dsd)/tan(alpha)-1;
            Tauv(i) = Tau;

        % Calculations of tooth deviations due to profile modification
            dx1 = x1-md1; % Travel distance from the contact point 1 on the line of action (useful) to the pinion (mm)
            dx2 = x2d-md2; % Travel distance from contact point 2 on the line of action (useful) to the crown (mm)

        % Stiffness from energy method
            uf1 = l1-h1*tan(phi1); Sf1 = thetab1*di1;
            uf2 = l2-h2*tan(phi2); Sf2 = thetab2*di2;
            uf1d = l1d-h1d*tan(phi1d); uf2d = l2d-h2d*tan(phi2d);
            Kf1 = (cos(phi1)^2/(E*b/1000))*(XL1*(uf1/Sf1)^2 + XM1*(uf1/Sf1) + XP1*(1+XQ1*tan(phi1)^2));
            Kf1d = (cos(phi1d)^2/(E*b/1000))*(XL1*(uf1d/Sf1)^2 + XM1*(uf1d/Sf1) + XP1*(1+XQ1*tan(phi1d)^2));
            Kf2 = (cos(phi2)^2/(E*b/1000))*(XL2*(uf2/Sf2)^2 + XM2*(uf2/Sf2) + XP2*(1+XQ2*tan(phi2)^2));
            Kf2d = (cos(phi2d)^2/(E*b/1000))*(XL2*(uf2d/Sf2)^2 + XM2*(uf2d/Sf2) + XP2*(1+XQ2*tan(phi2d)^2));

            if theta1 >= tdouble(1) && theta1 < tdouble(2) % At phase 1 (duplo)
                K1 = Ikbi(phi1,thetab1)+Ikai(phi1,thetab1)+Iksi(phi1,thetab1)+Kf1;
                K2 = Ikbi(phi2,thetab2)+Ikai(phi2,thetab2)+Iksi(phi2,thetab2)+Kf2;
                K1d = Ikbi(phi1d,thetab1)+Ikai(phi1d,thetab1)+Iksi(phi1d,thetab1)+Kf1d;
                K2d = Ikbi(phi2d,thetab2)+Ikai(phi2d,thetab2)+Iksi(phi2d,thetab2)+Kf2d;
                Kte1 = 1/(1/Kh+K1+K2); Kte2 = 1/(1/Kh+K1d+K2d);
                Kte(i) = Kte1 + Kte2;
                K12(i,1) = Kte1; K12(i,2) = Kte2;

                deltax1 = (deltamax/1000)*(double-dx1)/double; % Profile error due to the modification (with relation to gear 2 pair 1) (mm)
                deltax2 = deltax1*(double-dx2)/(double-dx1); % Profile error due to the modification (with relation to the gear 1 pair 2) (mm)
                delta(i,1) = deltax1/1000; delta(i,2) = deltax2/1000; % Error of the meshing pairs (mm)

                Ktem(i) = F*Kte(i)/(F+Kte1*deltax1/1000+Kte2*deltax2/1000); % Equivalent stiffness (N/m)

                R1(i) = Kte1/Kte(i);
            end

            if theta1 >= tsingle(1) && theta1 <= tsingle(2) % Esta em fase 2 (simples)
                K1 = Ikbi(phi1,thetab1)+Ikai(phi1,thetab1)+Iksi(phi1,thetab1)+Kf1;
                K2 = Ikbi(phi2,thetab2)+Ikai(phi2,thetab2)+Iksi(phi2,thetab2)+Kf2;
                Kte1 = 1/(1/Kh+K1+K2); Kte2 = 0;
                Kte(i) = Kte1 + Kte2;
                K12(i,1) = Kte1; K12(i,2) = Kte2;

                delta(i,1) = 0; delta(i,2) = 0;

                Ktem(i) = Kte(i); % Equivalent stiffness  (N/m)

                R1(i) = Kte1/Kte(i);
                mem = i;
            end

            if theta1 > (thetad+thetas) && theta1 <= (2*thetad+thetas) % Pahse 3 (double)
                R1(i) = 1- R1(i-mem);
            end

         % Meshing errors
            % Energy
                TE(i) = F/(Kte1+Kte2); % Theorical transmission error (m)
                dTE1 = (Kte1*deltax1/1000)/Kte(i); % Equivalent static transmission error (pair 1) (m)
                dTE2 = (Kte2*deltax2/1000)/Kte(i); % Equivalent static transmission error (pair 2) (m)
                dTE(i) = dTE1 + dTE2; % Total theorical static transmission error (m)
                TEm(i) = TE(i) + dTE(i); % Static transmission error (only caused by TPM) (m)

        % Damping
            Fn = F*cos(alphax1dsd); % Normal force at the tooth surface at the contact point (N)
            Ft = F*sin(alphax1dsd); % Tangential force at the tooth surface at the contact point (N)
            vrel = abs(w*x1dsd/1000-(rp1/rp2)*w*x2dsd/1000); vrelv(i) = vrel;
            mu = 0.1211; % 0.12*(Ft*Ra/(v0*vrel*rx1dsd/1000))^0.25; % Friction coefficient. Mean value mu = 0.015 / mu máx (primitive point) = 0.1211.
            Fa = Fn*mu; % Friction forc
            B0 = Fn/Kte(i); % Vibration amplitude (equivalent to the static deformation);
            Wa = 4*Fa*B0; % Work done by friction
            U = (Fn^2)/(2*Kte(i)); % Energy loss due friction (equivalent to the energy loss in the tooth deformation);
            zeta = Wa/(4*pi*U); % Damping factor

            C(i) = 2*zeta*sqrt(me*Kte(i)); % Damping coefficient (kg/s)
            Cm(i) = 2*zeta*sqrt(memod*Ktem(i)); % Equivalent damping coefficient (kg/s)
            
        i = i + 1; % Counter to organize the vectors

    end

end


