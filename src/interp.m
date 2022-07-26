function [tint,Tauint,Kteint,K1int,K2int,Ktmint,Cint,Cmint,Fint,d1int,d2int,KexData,KeyData,CxData,CyData] = interp(tempo,step,w,Tauv,Kte,K12,Ktem,C,Cm,R1,F,delta)

% TEMPO
        time = linspace(0,tempo,step*round(tempo*w/(2*pi)));
        num = linspace(1,step*round(tempo*w/(2*pi)),step*round(tempo*w/(2*pi)));
        [txData,tyData] = prepareCurveData(num,time);
        [tint,~] = fit(txData,tyData,'pchipinterp');

        Tauexp = repmat(Tauv,1,round(tempo*w/(2*pi)));
        [TauxData,TauyData] = prepareCurveData(Tauexp,time);
        [Tauint,~] = fit(TauxData,TauyData,'pchipinterp');

    % RIGIDEZ
        % Energia
            Kteexp = repmat(Kte,1,round(tempo*w/(2*pi)));
            [KexData,KeyData] = prepareCurveData(time,Kteexp);
            [Kteint,~] = fit(KexData,KeyData,'pchipinterp');

        % Separada
            K1exp = repmat(K12(:,1)',1,round(tempo*w/(2*pi)));
            [K1xData,K1yData] = prepareCurveData(time,K1exp);
            [K1int,~] = fit(K1xData,K1yData,'pchipinterp');
            
            K2exp = repmat(K12(:,2)',1,round(tempo*w/(2*pi)));
            [K2xData,K2yData] = prepareCurveData(time,K2exp);
            [K2int,~] = fit(K2xData,K2yData,'pchipinterp');

        % Modificada
            Ktmxp = repmat(Ktem,1,round(tempo*w/(2*pi)));
            [KmxData,KmyData] = prepareCurveData(time,Ktmxp);
            [Ktmint,~] = fit(KmxData,KmyData,'pchipinterp');

    % AMORTECIMENTO
        % Energia
            Cexp = repmat(C,1,round(tempo*w/(2*pi)));
            [CxData,CyData] = prepareCurveData(time,Cexp);
            [Cint,~] = fit(CxData,CyData,'pchipinterp');

        % Modificado
            Cmexp = repmat(Cm,1,round(tempo*w/(2*pi)));
            [CmxData,CmyData] = prepareCurveData(time,Cmexp);
            [Cmint,~] = fit(CmxData,CmyData,'pchipinterp');

    % FORCA
        Fexp = repmat(R1*F,1,round(tempo*w/(2*pi)));
        [FxData,FyData] = prepareCurveData(time,Fexp);
        [Fint,~] = fit(FxData,FyData,'pchipinterp');

    % QUANDITADE DE MODIFICACAO
        d1exp = repmat(delta(:,1)',1,round(tempo*w/(2*pi)));
        [d1xData,d1yData] = prepareCurveData(time,d1exp);
        [d1int,~] = fit(d1xData,d1yData,'pchipinterp');

        d2exp = repmat(delta(:,2)',1,round(tempo*w/(2*pi)));
        [d2xData,d2yData] = prepareCurveData(time,d2exp);
        [d2int,~] = fit(d2xData,d2yData,'pchipinterp');
end