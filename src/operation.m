function [w,w1,w2,T,F,freq,ciclo] = operation(Pot,W,db,dp1,dp2)
    w = W*(2*pi/60); % Pinion rotation (rad/s)
    w1 = w; w2 = -w1*dp1/dp2; % (rad/s)
    freq = w/(2*pi); % Frequency of pinion rotation (Hz)
    ciclo = 2*pi/w; % Period of revolution of the pinion (s) (used as the minimum analysis time)
    T = Pot*(60/(2*pi))/W; % Torque acting in the line of action (Nm) (pinion)
    F = T/(db/2000); % Force acting in the line of action (N)
end
