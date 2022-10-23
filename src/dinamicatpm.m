function dxdt = dinamicatpm(x,t,K,C,F,me,d1,d2,K1,K2) % PADRAO

%% MATRIZ DO SISTEMA E ENTRADA
A = [0 1; -K(t)/me -C(t)/me];
B = [0; 1/me];

%% SISTEMA
    dxdt = A*x + B*(F + K1(t)*d1(t) + K2(t)*d2(t));
end
