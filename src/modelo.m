function [t1,t2,estados,estadostpm,acel,aceltpm,tvmf,tvmftpm] = modelo(step,ciclos,tempo,Kteint,Ktmint,K1int,K2int,Cint,Cmint,d1int,d2int,F,me,memod,x0)

tspan = 0:(tempo/(ciclos*step)):tempo-(tempo/(ciclos*step)); %[0 tempo]; % Tempo de analise
opts = odeset('Refine',5,'Stats','on');

% SIMULACAO DO SISTEMA ORIGINAL
    [t1, estados] = ode15s(@(t,x) dinamica(x,t,Kteint,Cint,F,me),tspan,x0);

    % Calculo da aceleracao
        acel = zeros(size(estados)); tvmf = zeros(size(estados)); clear k;
        for k = 1:size(estados,1)
            acel(k,:) = aceleracao(estados(k,:),Kteint(t1(k)),Cint(t1(k)),F,me);
            tvmf(k,:) = Kteint(t1(k))*estados(k,1) + Cint(t1(k))*estados(k,2);% + me*acel(k,2);
        end
        acel = acel(:,2); tvmf = tvmf(:,1);

% SIMULACAO DO SISTEMA COM TPM
     [t2, estadostpm] = ode15s(@(t,x) dinamicatpm(x,t,Kteint,Cint,F,memod,d1int,d2int,K1int,K2int),tspan,x0);

    % Calculo da aceleracao
        aceltpm = zeros(size(estadostpm)); tvmftpm = zeros(size(estadostpm)); clear k;
        for k = 1:size(estadostpm,1)
            aceltpm(k,:) = aceleracaotpm(estadostpm(k,:),Kteint(t2(k)),Cint(t2(k)),F,memod,d1int(t2(k)),d2int(t2(k)),K1int(t2(k)),K2int(t2(k)));
            tvmftpm(k,:) = Kteint(t2(k))*estadostpm(k,1) + Cint(t2(k))*estadostpm(k,2) - K1int(t2(k))*d1int(t2(k)) - K2int(t2(k))*d2int(t2(k)); % + memod*aceltpm(k,2);
        end
        aceltpm = aceltpm(:,2); tvmftpm = tvmftpm(:,1);
end
