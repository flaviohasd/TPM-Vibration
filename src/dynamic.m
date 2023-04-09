function dxdt = dynamic(x,t,K,C,F,me)

%% SYSTEM MATRIX AND INPUT MATRIX
A = [0 1; -K(t)/me -C(t)/me];
B = [0; 1/me];

%% SYSTEM
    dxdt = A*x + B*F;
end