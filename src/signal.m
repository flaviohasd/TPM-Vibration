function [RMS,RMStpm,Amp,Amptpm,Espectro,Espectrotpm,ft,fttpm] = signal(states,statestpm,freq)
    step = size(states(:));

    RMS = rms(states(step/2:step)); % Mean quadratic value of the original system
    RMStpm = rms(statestpm(step/2:step)); % VMean quadratic value of the modified system

    Amp = peak2peak(states(step/2:step)); % Peak-to-peak amplitude of the original system
    Amptpm = peak2peak(statestpm(step/2:step)); % Peak-to-peak amplitude of the modified system

    Espectro = envspectrum(states(step/2:step),freq); % Frequency spectrum of the original signal
    Espectrotpm = envspectrum(statestpm(step/2:step),freq); % Frequency spectrum of the modified signal
    
    pot2 = 2^nextpow2(length(states(:,1)));
    ft = fft(states(:,1),pot2);
    fttpm = fft(statestpm(:,1),pot2);
 end

