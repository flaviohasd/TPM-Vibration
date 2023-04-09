function [RMS,RMStpm,Amp,Amptpm,Espectro,Espectrotpm,ft,fttpm] = signal(estados,estadostpm,freq)

    RMS = rms(estados(size(estados(:,1),1)/2:size(estados(:,1),1),1)); % Mean quadratic value of the original system
    RMStpm = rms(estadostpm(size(estadostpm(:,1),1)/2:size(estadostpm(:,1),1),1)); % VMean quadratic value of the modified system

    Amp = peak2peak(abs(estados(size(estados(:,1),1)/2:size(estados(:,1),1),1))); % Peak-to-peak amplitude of the original system
    Amptpm = peak2peak(abs(estadostpm(size(estadostpm(:,1),1)/2:size(estadostpm(:,1),1),1))); % Peak-to-peak amplitude of the modified system

    Espectro = envspectrum(estados(size(estados(:,1),1)/2:size(estados(:,1),1),1),freq); % Frequency spectrum of the original signal
    Espectrotpm = envspectrum(estadostpm(size(estadostpm(:,1),1)/2:size(estadostpm(:,1),1),1),freq); % Frequency spectrum of themodified signal
    
    pot2 = 2^nextpow2(length(estados(:,1)));
    ft = fft(estados(:,1),pot2);
    fttpm = fft(estadostpm(:,1),pot2);
 end

