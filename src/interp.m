function [tint,Kteint,K1int,K2int,Ktmint,Cint,Cmint,Fint,d1int,d2int,KexData,KeyData,CxData,CyData,TEexp,TEmexp] = interp(ciclos,tempo,step,w,Tauv,Kte,K12,Ktem,C,Cm,R1,F,delta,TE,TEm)

% TEMPO
        time = linspace(0,tempo,step*ciclos);
        num = linspace(1,step*ciclos,step*ciclos);
        [txData,tyData] = prepareCurveData(num,time);
        [tint,~] = fit(txData,tyData,'pchipinterp');

    % RIGIDEZ
        % Energia
            Kteexp = repmat(Kte,1,ciclos);
            [KexData,KeyData] = prepareCurveData(time,Kteexp);
            [Kteint,~] = fit(KexData,KeyData,'pchipinterp');

        % Separada
            K1exp = repmat(K12(:,1)',1,ciclos);
            [K1xData,K1yData] = prepareCurveData(time,K1exp);
            [K1int,~] = fit(K1xData,K1yData,'pchipinterp');
            
            K2exp = repmat(K12(:,2)',1,ciclos);
            [K2xData,K2yData] = prepareCurveData(time,K2exp);
            [K2int,~] = fit(K2xData,K2yData,'pchipinterp');

        % Modificada
            Ktmxp = repmat(Ktem,1,ciclos);
            [KmxData,KmyData] = prepareCurveData(time,Ktmxp);
            [Ktmint,~] = fit(KmxData,KmyData,'pchipinterp');

    % AMORTECIMENTO
        % Energia
            Cexp = repmat(C,1,ciclos);
            [CxData,CyData] = prepareCurveData(time,Cexp);
            [Cint,~] = fit(CxData,CyData,'pchipinterp');

        % Modificado
            Cmexp = repmat(Cm,1,ciclos);
            [CmxData,CmyData] = prepareCurveData(time,Cmexp);
            [Cmint,~] = fit(CmxData,CmyData,'pchipinterp');

    % FORCA
        Fexp = repmat(R1*F,1,ciclos);
        [FxData,FyData] = prepareCurveData(time,Fexp);
        [Fint,~] = fit(FxData,FyData,'pchipinterp');

    % QUANDITADE DE MODIFICACAO
        d1exp = repmat(delta(:,1)',1,ciclos);
        [d1xData,d1yData] = prepareCurveData(time,d1exp);
        [d1int,~] = fit(d1xData,d1yData,'pchipinterp');

        d2exp = repmat(delta(:,2)',1,ciclos);
        [d2xData,d2yData] = prepareCurveData(time,d2exp);
        [d2int,~] = fit(d2xData,d2yData,'pchipinterp');
        
    % ERROS
        TEexp = repmat(TE,1,ciclos);
        TEmexp = repmat(TEm,1,ciclos);
        
end