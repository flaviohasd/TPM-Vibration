%% OTIMIZACAO
options = optimset('PlotFcns',@optimplotfval,'Display','iter','TolX',1e-4);
tic

deltamaxopt = fminbnd(@simulacao,0,100,options);

disp(['Tempo para realizar a otimização: ',num2str(toc/60),' minutos'])
disp(' ')
disp(['Valor de deltamax que minimiza a vibração: ',num2str(deltamaxopt)])

mm = simulacao2(deltamaxopt);

% Luo 38.59 (2k - 80k)
% Tian 93.7337 (2k - 80k)