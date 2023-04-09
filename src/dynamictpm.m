function dxdt = dynamictpm(x,t,K,C,F,me,d1,d2,K1,K2)

%% SYSTEM MATRIX AND INPUT MATRIX
A = [0 1; -K(t)/me -C(t)/me];
B = [0; 1/me];

%% SYSTEM
    dxdt = A*x + B*(F + K1(t)*d1(t) + K2(t)*d2(t));
end