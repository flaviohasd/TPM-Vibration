function dxdt = dinamica(x,t,K,C,F,me) % PADRAO

%% MATRIZ DO SISTEMA E ENTRADA
A = [0 1; -K(t)/me -C(t)/me];
B = [0; 1/me];

%% SISTEMA
    dxdt = A*x + B*F;
end
