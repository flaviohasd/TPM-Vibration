function [RMS,RMStpm,Amp,Amptpm,Espectro,Espectrotpm,ft,fttpm] = sinal(estados,estadostpm,freq)

    RMS = rms(estados(size(estados(:,1),1)/2:size(estados(:,1),1),1)); % Valor quadratico medio do sistema original
    RMStpm = rms(estadostpm(size(estadostpm(:,1),1)/2:size(estadostpm(:,1),1),1)); % Valor quadratico medio do sistema modificado

    Amp = max(abs(estados(size(estados(:,1),1)/2:size(estados(:,1),1),1))); % Amplitude pico do sistema original
    Amptpm = max(abs(estadostpm(size(estadostpm(:,1),1)/2:size(estadostpm(:,1),1),1))); % Amplitude pico do sistema modificado

    Espectro = envspectrum(estados(size(estados(:,1),1)/2:size(estados(:,1),1),1),freq); % Espectro da frequencia do sinal original
    Espectrotpm = envspectrum(estadostpm(size(estadostpm(:,1),1)/2:size(estadostpm(:,1),1),1),freq); % Espectro da frequencia do sinal modificado
    
    pot2 = 2^nextpow2(length(estados(:,1)));
    ft = fft(estados(:,1),pot2);
    fttpm = fft(estadostpm(:,1),pot2);
 end

