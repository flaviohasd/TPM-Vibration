function [t1,t2,states,statestpm,acel,aceltpm,dmf,dmftpm] = model(step,cycles,time,Kteint,Ktmint,K1int,K2int,Cint,Cmint,d1int,d2int,F,me,memod,x0)

tspan = 0:(time/(cycles*step)):time-(time/(cycles*step)); %[0 time]; % Time of analysis
opts = odeset('Refine',5,'Stats','on');

% SIMULATION OF THE ORIGINAL SYSTEM
    [t1, states] = ode15s(@(t,x) dynamic(x,t,Kteint,Cint,F,me),tspan,x0);

    % Acceleration calculation
        acel = zeros(size(states)); dmf = zeros(size(states)); clear k;
        for k = 1:size(states,1)
            acel(k,:) = accel(states(k,:),Kteint(t1(k)),Cint(t1(k)),F,me);
            dmf(k,:) = Kteint(t1(k))*states(k,1) + Cint(t1(k))*states(k,2);
        end
        acel = acel(:,2); dmf = dmf(:,1);

% SIMULATION OF THE MODIFIED SYSTEM
     [t2, statestpm] = ode15s(@(t,x) dynamictpm(x,t,Kteint,Cint,F,memod,d1int,d2int,K1int,K2int),tspan,x0);

    % Acceleration and force calculation
        aceltpm = zeros(size(statestpm)); dmftpm = zeros(size(statestpm)); clear k;
        for k = 1:size(statestpm,1)
            aceltpm(k,:) = acceltpm(statestpm(k,:),Kteint(t2(k)),Cint(t2(k)),F,memod,d1int(t2(k)),d2int(t2(k)),K1int(t2(k)),K2int(t2(k)));
            dmftpm(k,:) = Kteint(t2(k))*statestpm(k,1) + Cint(t2(k))*statestpm(k,2) - K1int(t2(k))*d1int(t2(k)) - K2int(t2(k))*d2int(t2(k));
        end
        aceltpm = aceltpm(:,2); dmftpm = dmftpm(:,1);
end
