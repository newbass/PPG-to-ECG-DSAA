function [PPG1, PPG2, ACCX, ACCY, ACCZ, PPGave, srate, Freq_Range, PPGaveFFT, ACCXFFT, ACCYFFT, ACCZFFT] = data_reconfig(PPG1, PPG2, ACCX, ACCY, ACCZ, b_coeff, a_coeff, FFTres, Low_freq, High_freq)
    % filtering
    PPGave = 0.5*(PPG1-mean(PPG1))/(std(PPG1))+0.5*(PPG2-mean(PPG2))/(std(PPG2)); % mean overall
    [ACCX, ACCY, ACCZ, PPG1, PPG2] = deal(filter(b_coeff,a_coeff,ACCX), filter(b_coeff,a_coeff,ACCY), filter(b_coeff,a_coeff,ACCZ), filter(b_coeff,a_coeff,PPG1), filter(b_coeff,a_coeff,PPG2));

    % downsampling to 25Hz
    [ACCX, ACCY, ACCZ, PPGave] = deal(downsample(ACCX,5), downsample(ACCY,5), downsample(ACCZ,5), downsample(PPGave,5));
    srate = 25; % new sampling rate

    % Periodogram
    PPGaveFFT = fft(PPGave,FFTres);
    Freq_Range = linspace(0,srate,size(PPGaveFFT,2));

    % finding the indices for the range of interest
    [~,lowR] = (min(abs(Freq_Range - High_freq)));
    [~,highR] = (min(abs(Freq_Range - Low_freq)));

    %  Getting rid of most spectra outside the range of interest
    [ACCXFFT, ACCYFFT, ACCZFFT] = deal(fft(ACCX,FFTres), fft(ACCY,FFTres), fft(ACCZ,FFTres));
    [ACCXFFT, ACCYFFT, ACCZFFT, Freq_Range, PPGaveFFT] = deal(ACCXFFT(lowR:highR), ACCYFFT(lowR:highR), ACCZFFT(lowR:highR), Freq_Range(lowR:highR), PPGaveFFT(lowR:highR));
end