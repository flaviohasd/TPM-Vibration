function [t,states,statestpm,acel,aceltpm,dmf,dmftpm] = model(step,cycles,time,Kteint,Ktmint,K1int,K2int,Cint,Cmint,d1int,d2int,F,me,memod,x0,w)

tspan = 0:(time/(cycles*step)):time-(time/(cycles*step)); %[0 time]; % Time of analysis
t = linspace(0,(2*pi)/w,step); % Step time
opts = odeset('Refine',5,'Stats','on');

% SIMULATION OF THE ORIGINAL SYSTEM
    [~, states] = ode15s(@(t,x) dynamic(x,t,Kteint,Cint,F,me),tspan,x0);
 
    % Acceleration calculation
        acel = zeros(size(states)); dmf = zeros(size(states)); clear k;
        for k = 1:size(states,1)
            acel(k,:) = accel(states(k,:),Kteint(t(k)),Cint(t(k)),F,me);
            dmf(k,:) = Kteint(t(k))*states(k,1) + Cint(t(k))*states(k,2);
        end
        acel = acel(:,2); dmf = dmf(:,1);

% SIMULATION OF THE MODIFIED SYSTEM
     [~,statestpm] = ode15s(@(t,x) dynamictpm(x,t,Kteint,Cint,F,memod,d1int,d2int,K1int,K2int),tspan,x0);

    % Acceleration and force calculation
        aceltpm = zeros(size(statestpm)); dmftpm = zeros(size(statestpm)); clear k;
        for k = 1:size(statestpm,1)
            aceltpm(k,:) = acceltpm(statestpm(k,:),Kteint(t(k)),Cint(t(k)),F,memod,d1int(t(k)),d2int(t(k)),K1int(t(k)),K2int(t(k)));
            dmftpm(k,:) = Kteint(t(k))*statestpm(k,1) + Cint(t(k))*statestpm(k,2) - K1int(t(k))*d1int(t(k)) - K2int(t(k))*d2int(t(k));
        end
        aceltpm = aceltpm(:,2); dmftpm = dmftpm(:,1);
end
