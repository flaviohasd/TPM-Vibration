%% OTIMIZACAO
options = optimset('PlotFcns',@optimplotfval,'Display','iter','TolX',1e-5);
tic
deltamaxopt = fminbnd(@simulacao,0,50,options);
disp(['Valor de deltamax que minimiza a vibração: ',num2str(deltamaxopt)])
toc