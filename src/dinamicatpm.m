function dxdt = dinamicatpm(x,t,K,C,F,me,d1,d2,K1,K2) % PADRAO

%% MATRIZ DO SISTEMA E ENTRADA
A = [0 1; -K(t)/me -C(t)/me];
B = [0; 1/me];

%% SISTEMA
    dxdt = A*x + B*(F + K1(t)*d1(t) + K2(t)*d2(t));
end

%         [t, estados] = ode15s(@(t,x) dinamica(x,t,Kteint,Cint,F,rb1/1000,rb2/1000,J1,J2),tspan,x0,opts); % Expandido (Mesmo resultado da SS MATRIZ)
%         [t, estados] = ode15s(@(t) dinamica(time,Kteint,Cint,F,rb1,rb2,J1,J2),tspan,x0,opts); % Modal (Ruim)
%         [t, estados] = ode15s(@(t,x) dinamica(x,t,Kteint,Cint,[rb1*F/1000;-rb2*F/1000],inv([J1,0;0,J2]),rb1/1000,rb2/1000),tspan,x0,opts); % SS MATRIZ (Estranho)



% function dxdt = dinamica(x,t,tempo,Ktexp,Cexp,F,me) % COM INTERPOLACAO INTERNA
% 
% %% INTERPOLACAO DAS VARIAVEIS NO TEMPO
%     % RIGIDEZ
%         [KxData,KyData] = prepareCurveData(tempo,Ktexp);
%         [Ktint,~] = fit(KxData,KyData,'pchipinterp');
% 
%     % AMORTECIMENTO
%         [CxData,CyData] = prepareCurveData(tempo,Cexp);
%         [Cint,~] = fit(CxData,CyData,'pchipinterp');
% 
% %% MATRIZ DO SISTEMA E ENTRADA
% A = [0 1; -Ktint(t)/me -Cint(t)/me];
% B = [0;1];
% 
% %% SISTEMA
%     dxdt = A*x+B*F;
% end



% function dxdt = dinamica(x,t,Ktint,Cint,F,rb1,rb2,J1,J2) % EXPANDIDO
% %% MATRIZ DO SISTEMA E ENTRADA
% A = [           0,                   1,                  0,                    0;
%      -(Ktint(t)*rb1^2/J1), -(Cint(t)*rb1^2/J1), (Ktint(t)*rb1*rb2/J1), (Cint(t)*rb1*rb2);
%                 0,                   0,                  0,                  1;
%       (Ktint(t)*rb1*rb2/J2), (Cint(t)*rb1*rb2/J2), -(Ktint(t)*rb2^2/J2), -(Cint(t)*rb2^2/J2)];
% 
% B = [0; (rb1/J1); 0; -(rb2/J2)];
% 
% %% SISTEMA
%     dxdt = A*x+B*F;
% end



% function dxdt = dinamica(x,t,Ktint,Cint,F,iM,rb1,rb2) % SS MATRIZ
% 
% %% MATRIZ DO SISTEMA E ENTRADA
% K = [Ktint(t)*rb1^2,-Ktint(t)*rb1*rb2;-Ktint(t)*rb1*rb2,Ktint(t)*rb2^2];
% C = [Cint(t)*rb1^2,-Cint(t)*rb1*rb2;-Cint(t)*rb1*rb2,Cint(t)*rb2^2];
% 
% A = [0 0 1 0; 0 0 0 1; -iM*K -iM*C];
% B = [0; 0; -iM*F];
% 
% %% SISTEMA
%     dxdt = A*x + B;
% end



% function xva = dinamica(Ktint,Cint,F,rb1,rb2,J1,J2) % MODAL
%     syms t
% 
%         % Modelo dinamico
%         M = [J1, 0;
%             0, J2];
%         
%         F = [rb1*F;
%              -rb2*F];
%          
%          Z = [0.1;
%               0.1];     
%           
%          L = chol(M); Li = inv(L);
%          
%         i = 1;
%         for te = 0:0.001:tempo
%             syms t
% 
%             C = [Cint(te)*rb1^2, -Cint(te)*rb1*rb2;
%                  -Cint(te)*rb1*rb2, Cint(te)*rb2^2];
% 
%             K = [Ktint(te)*rb1^2, -Ktint(te)*rb1*rb2;
%                  -Ktint(te)*rb1*rb2, Ktint(te)*rb2^2];
% 
%         % ANALISE MODAL
% 
%             Ktil = Li*K*Li; % Li*K*Li
%             Ctil = Li*C*Li; % Li*C*Li
%             [P, Delta] = eig(Ktil);
%             Pt = transpose(P);
%             S = Li*P;
%             Si = Pt*L;
%             Cr = Pt*Ctil*P;
%             Kr = Delta;
%             Fr = Pt*Li*F; % Pt*Li*F
%             W = sqrt(Delta); Wd = W.*sqrt(1.-Z.^2);
% 
%         % Resposta Modal
% 
%             r1 = Fr(1)*t; % Frequencia natual w1 = 0 -> Resposta por laplace
%             Amp = (Fr(2)/Wd(2,2));
%             Exp = -diag(Z(2))*W(2,2).*t;
%             Rm = [r1; Amp*exp(Exp)*sin(Wd(2,2)*t)]; % Para condicoes iniciais nulas
% 
%             % Sistema Fisico
%             X = vpa(Li*P*Rm,4);
%             x = vpa(X,4); v = diff(x,t); a = diff(v,t);
%             x = subs(x,t,te); v = subs(v,t,te); a = subs(a,t,te);
% 
%             p1(i) = x(1); p2(i) = x(2); v1(i) = v(1); v2(i) = v(2); a1(i) = a(1); a2(i) = a(2);
%             
%             i = i + 1;
%         end
% end