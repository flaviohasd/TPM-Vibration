function [tint,Kteint,K1int,K2int,Ktmint,Cint,Cmint,Fint,d1int,d2int,KexData,KeyData,CxData,CyData,TEexp,TEmexp] = interpolation(cycles,time,step,w,Tauv,Kte,K12,Ktem,C,Cm,R1,F,delta,TE,TEm)

% TIME
        time = linspace(0,time,step*cycles);
        num = linspace(1,step*cycles,step*cycles);
        [txData,tyData] = prepareCurveData(num,time);
        [tint,~] = fit(txData,tyData,'pchipinterp');

    % STIFFNESS
        % Energy
            Kteexp = repmat(Kte,1,cycles);
            [KexData,KeyData] = prepareCurveData(time,Kteexp);
            [Kteint,~] = fit(KexData,KeyData,'pchipinterp');

        % Separated
            K1exp = repmat(K12(:,1)',1,cycles);
            [K1xData,K1yData] = prepareCurveData(time,K1exp);
            [K1int,~] = fit(K1xData,K1yData,'pchipinterp');
            
            K2exp = repmat(K12(:,2)',1,cycles);
            [K2xData,K2yData] = prepareCurveData(time,K2exp);
            [K2int,~] = fit(K2xData,K2yData,'pchipinterp');

        % Modified
            Ktmxp = repmat(Ktem,1,cycles);
            [KmxData,KmyData] = prepareCurveData(time,Ktmxp);
            [Ktmint,~] = fit(KmxData,KmyData,'pchipinterp');

    % DAMPING
        % Energy
            Cexp = repmat(C,1,cycles);
            [CxData,CyData] = prepareCurveData(time,Cexp);
            [Cint,~] = fit(CxData,CyData,'pchipinterp');

        % Modified
            Cmexp = repmat(Cm,1,cycles);
            [CmxData,CmyData] = prepareCurveData(time,Cmexp);
            [Cmint,~] = fit(CmxData,CmyData,'pchipinterp');

    % FORCE
        Fexp = repmat(R1*F,1,cycles);
        [FxData,FyData] = prepareCurveData(time,Fexp);
        [Fint,~] = fit(FxData,FyData,'pchipinterp');

    % QUALITY OF INTERPOLATION
        d1exp = repmat(delta(:,1)',1,cycles);
        [d1xData,d1yData] = prepareCurveData(time,d1exp);
        [d1int,~] = fit(d1xData,d1yData,'pchipinterp');

        d2exp = repmat(delta(:,2)',1,cycles);
        [d2xData,d2yData] = prepareCurveData(time,d2exp);
        [d2int,~] = fit(d2xData,d2yData,'pchipinterp');
        
    % ERRORS
        TEexp = repmat(TE,1,cycles);
        TEmexp = repmat(TEm,1,cycles);
        
end