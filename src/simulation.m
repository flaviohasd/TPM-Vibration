function RMSatpm = simulation(deltamax) % Perform the optimization using the RMS of the acceleration signal
    syms rx1 rx2 alpha1 alpha2 x

  %% PARAMETERS OF THE CODE
        cycle = 1; % Cycle (gear revolution) quantity for analysis (min=1, bigger it is -> slower the code and it may cause interpolation errors)
        step = 1e4; % Quantity of iteractions to use in the loop for the stiffness calculation. (OBS: Use at least 1e4 to avoid big iteraction errors)
        % deltamax = 43.5334; % Maximum amount of modification (µm)      %% Leave it commented for the optimization process.
        x0 = [0; 0]; % Initial conditions of the space-state model: [displacement(0) velocity(0)] of the dynamic model

    %% CONSTANTS
        % GEOMETRY AND PROPERTY OF GEARS
            z1 = 27; % N pinhon teeth
            z2 = 35; % N gear teerh
            m = 3; % Normal modulus (mm) (=25.4/P or P = [1/in]) 
            alpha = deg2rad(20); % Pressure angle (rad)
            b = 25; % Width (mm)

            E = 206e9; % Elasticity Modulus (Pa)
            v = 0.3; % Poisson's coefficient
            p = 7850; % Material density of gears (kg/m^3)
            Ra = 0.63/1e6; % Surface roughtness (m)

        % OPERATING CONDITIONS
            Pow = 80e3; % Input power (W) (pinion)
            W = 2000; % Input rotation (rpm) (pinion)

            p0 = 870;  % Lubricant density (kg/m^3)
            t0 = 60; % Operating temperature (°C)
            v0 = 38.5+(t0-100)*(38.5-320)/(100-40); % Kinematic viscosity (cst) -  320 (40°C) / 38.5 (100°C)
            
            % Inconsistency check
                % Pinion gear
                    if z1 > z2
                        temp = z1;
                        z1 = z2;
                        z2 = temp;
                    end

                % Pressure angle
                    alphacrit = asin(sqrt(2/z1));   % Testing the critical pressure angle
                    if alphacrit > alpha
                        fprintf('Pressure angle less than the critial angle: %.2f°! - Solution: increase the number of teeth. \n',rad2deg(alphacritdeg))
                        return
                    end
                    
                % Number of teeth
                    zmin = ceil((2/(3*(sin(alpha))^2))*(1+sqrt(1+3*sin(alpha)^2)));
                    if z1 < zmin
                        fprintf('Number of teeth is less than minimum: %.2f°! - Solution: increase the number of teeth. \n',zmin)
                    return
                    end

        % DIMENTIONS
            [dp1,de1,db1,di1,rd1,ra1,rb1,ri1,rp1,rf1,hf1,Sa1,Sd1,hr1,h1,hi1,dp2,de2,db2,di2,rd2,ra2,rb2,ri2,rp2,rf2,hf2,Sa2,Sd2,hr2,h2,hi2] = gears(z1,z2,m,alpha);

        % VELOCITY AND FORCE
            [w,w1,w2,T,F,freq,cyclo] = operation(Pow,W,db1,dp1,dp2); time = cyclo*cycle;

    % MESHING RELATIONS
       [pb,L,cr,PTH,md1,md2,thetad,thetas,alpha10,alpha20,thetab1,thetab2,double,single] = mesh(alpha,z1,z2,dp1,db1,ra1,rb1,dp2,ra2,rb2);

        if cr > 2 || cr < 1     %% Testing the contact ration criteria
            disp('Contact ratio outside of the valid interval!')
        end
    %% MASS AND MOMENT OF INERTIA
        [me,memod] = mass(z1,z2,p,deltamax,cr,pb,b,rb1,rb2,rf1,rf2,ri1,ri2,Sd1,Sd2,hr1,hr2,hi1,hi2,h1,h2,Sa1,Sa2);

    %% STIFFNESS AND DAMPING
        [Ikbi,Iksi,Ikai,Kh] = integr(E,b,v);

        [theta1vg,vrelv,R1,delta,Tauv,Kte,Ktem,K12,C,Cm,TE,TEm,dTE] = energy(step,deltamax,F,double,alpha,pb,z1,z2,di1,di2,rb1,rb2,rp1,rp2,w,E,b,thetad,thetas,thetab1,thetab2,alpha10,alpha20,L,hf1,hf2,md1,md2,Ikbi,Ikai,Iksi,Kh,me,memod,Ra,v0);

        % INTERPOLATION OF "K" AND "C" GIVEN THE TIME ENTRY
            [tint,Kteint,K1int,K2int,Ktmint,Cint,Cmint,Fint,d1int,d2int,KexData,KeyData,CxData,CyData,TEexp,TEmexp] = interpolation(cycle,time,step,w,Tauv,Kte,K12,Ktem,C,Cm,R1,F,delta,TE,TEm);

    %% DYNAMIC MODEL
        % SPACE-STATE SYSTEM
            [t1,t2,states,statestpm,acel,aceltpm,dmf,dmftpm] = model(step,cycle,time,Kteint,Ktmint,K1int,K2int,Cint,Cmint,d1int,d2int,F,me,memod,x0);

        % SIGNAL ANALYSIS
            % Steady-state
                % Displacement
                    [RMSx,RMSxtpm,Ampx,Ampxtpm,Espectrox,Espectroxtpm,ft,fttpm] = signal(states,statestpm,freq);
                    Fcx = Ampxtpm/RMSxtpm;

                % Acceleration
                    [RMSa,RMSatpm,Ampa,Ampatpm,Espectroa,Espectroatpm,ft,fttpm] = signal(acel,aceltpm,freq);

                % Force
                    [RMSf,RMSftpm,Ampf,Ampftpm,Espectrof,Espectroftpm,ft,fttpm] = signal(dmf,dmftpm,freq);

    %% PLOTS
         [mm] = plots(Ampx,Ampa,RMSx,RMSa,Ampxtpm,Ampatpm,RMSxtpm,RMSatpm,Kte,deltamax,w,theta1vg,C,vrelv,thetad,thetas,R1,Kteint,KexData,KeyData,Cint,CxData,CyData,t1,t2,states,statestpm,time,acel,aceltpm,dmf,dmftpm,Tauv,Ktem,TEexp,TEmexp,step,z1,delta,ft,fttpm);

end
