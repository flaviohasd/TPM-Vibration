function [w,w1,w2,T,F,freq,ciclo] = operacao(Pot,W,db,dp1,dp2)
    w = W*(2*pi/60); % Rotacao do pinhao (rad/s)
    w1 = w; w2 = -w1*dp1/dp2; % (rad/s)
    freq = w/(2*pi); % Frequencia de rotacao do pinhao (Hz)
    ciclo = 2*pi/w; % Periodo de revolucao do pinhao (s) (utilizado como tempo minimo de analise)
    T = Pot*(60/(2*pi))/W; % Torque atuante na linha de acao (Nm) (pinhao)
    F = T/(db/2000); % Forca atuante na linha de acao (N)
end
