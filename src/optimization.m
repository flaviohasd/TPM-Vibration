%% OPTIMIZATION
opt_tol = 1e-4; % Tolerance of the optimization process (algorithm criteria to stop)

options = optimset('PlotFcns',@optimplotfval,'Display','iter','TolX',opt_tol);
tic

deltamaxopt = fminbnd(@simulation,30,50,options);

disp(['Time to perform the optimizaiton: ',num2str(toc/60),' minutes'])
disp(' ')
disp(['Value of deltamax that minimizes the vibration: ',num2str(deltamaxopt)])

set(0,'defaultfigurecolor',[1 1 1]); mm = simulation2(deltamaxopt); 

load train; sound(y,Fs); %handel
