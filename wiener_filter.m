function [W1FFT, W11FFT, W2FFT, W21FFT, PPGaveFFTFIN] = wiener_filter(i, WFlength, PPGaveFFT, ACCXFFT, ACCYFFT, ACCZFFT, W1FFT, W11FFT, W2FFT, W21FFT, PPGaveFFTFIN)
    % Wiener filtering PPG-ACC, two types

    WC1 = WFlength; WC2 = WFlength;

    %Wiener 1, abs & normalised
    W1FFT(i,:) = (abs(PPGaveFFT))/max(abs(PPGaveFFT));
    if i==1
        W1PPGaveFFTALL = W1FFT(i,:); 
    else
        W1PPGaveFFTALL = mean(W1FFT(max(1,i-WC1):i,:),1);
    end
    W1PPGaveFFTALLnorm = (W1PPGaveFFTALL)/max(W1PPGaveFFTALL);
    [W1ACCXFFTnormal, W1ACCYFFTnormal, W1ACCZFFTnormal]= deal((abs(ACCXFFT))/max(abs(ACCXFFT)), (abs(ACCYFFT))/max(abs(ACCYFFT)), (abs(ACCZFFT))/max(abs(ACCZFFT)));
    WF1 = (1 - 1/3*(W1ACCXFFTnormal+W1ACCYFFTnormal+W1ACCZFFTnormal)./(W1PPGaveFFTALLnorm)); 
    WF1 (WF1<0) = -1; % limit negative -inf to -1
    W1PPGaveFFTClean(i,:) = abs(PPGaveFFT).*WF1;

    % Wiener 1, power2, works better on its own but not in ensemble
    W11FFT(i,:) = abs(PPGaveFFT).^2;
    if i==1
        W11PPGaveFFTALL = W11FFT(i,:); 
    else
        W11PPGaveFFTALL = mean(W11FFT(max(1,i-WC1):i,:),1); 
    end
    W11PPGaveFFTALLnorm = (W11PPGaveFFTALL)/max(W11PPGaveFFTALL);
    [W11ACCXFFTnormal, W11ACCYFFTnormal, W11ACCZFFTnormal] = deal((abs(ACCXFFT).^2)/max(abs(ACCXFFT.^2)), (abs(ACCYFFT).^2)/max(abs(ACCYFFT).^2), (abs(ACCZFFT).^2)/max(abs(ACCZFFT).^2));
    WF11 = (1 - 1/3*(W11ACCXFFTnormal+W11ACCYFFTnormal+W11ACCZFFTnormal)./(W11PPGaveFFTALLnorm)); 
    W11PPGaveFFTClean(i,:) = abs(PPGaveFFT).*(WF11);

    %Wiener 2, abs & normalised
    W2FFT(i,:) = (abs(PPGaveFFT))/max(abs(PPGaveFFT));
    if i==1
        W2PPGaveFFTALL = W2FFT(i,:); 
    else
        W2PPGaveFFTALL = mean(W2FFT(max(1,i-WC2):i,:),1); 
    end
    W2PPGaveFFTALLnormal = (W2PPGaveFFTALL)/max(W2PPGaveFFTALL);
    [W2ACCXFFTnormal, W2ACCYFFTnormal, W2ACCZFFTnormal] = deal((abs(ACCXFFT))/max(abs(ACCXFFT)), (abs(ACCYFFT))/max(abs(ACCYFFT)), (abs(ACCZFFT))/max(abs(ACCZFFT)));
    WF2 = W2PPGaveFFTALLnormal./(((W2ACCXFFTnormal+W2ACCYFFTnormal+W2ACCZFFTnormal)/3)+W2PPGaveFFTALLnormal); 
    W2PPGaveFFTClean(i,:) = abs(PPGaveFFT).*WF2;
    W2FFT(i,:) = (W2PPGaveFFTClean(i,:))/max(W2PPGaveFFTClean(i,:));

    %Wiener 2, power2, work better on its own, but not in ensemble
    W21FFT(i,:) = abs(PPGaveFFT).^2;
    if i==1
        W21PPGaveFFTALL = W21FFT(i,:); 
    else
        W21PPGaveFFTALL = mean(W21FFT(max(1,i-WC2):i,:),1); 
    end
    W21PPGaveFFTALLnormal = W21PPGaveFFTALL/max(W21PPGaveFFTALL);
    [W21ACCXFFTnormal, W21ACCYFFTnormal, W21ACCZFFTnormal] = deal(abs(ACCXFFT).^2/max(abs(ACCXFFT).^2), abs(ACCYFFT).^2/max(abs(ACCYFFT).^2), abs(ACCZFFT).^2/max(abs(ACCZFFT).^2));
    WF21 = W21PPGaveFFTALLnormal./(((W21ACCXFFTnormal+W21ACCYFFTnormal+W21ACCZFFTnormal)/3)+W21PPGaveFFTALLnormal);
    W21PPGaveFFTClean(i,:) = abs(PPGaveFFT).*WF21;
    W21FFT(i,:) = W21PPGaveFFTClean(i,:).^2;

    [W1PPGaveFFTClean(i,:), W11PPGaveFFTClean(i,:), W2PPGaveFFTClean(i,:), W21PPGaveFFTClean(i,:)] = deal(W1PPGaveFFTClean(i,:)/std(W1PPGaveFFTClean(i,:)), W11PPGaveFFTClean(i,:)/std(W11PPGaveFFTClean(i,:)), W2PPGaveFFTClean(i,:)/std(W2PPGaveFFTClean(i,:)), W21PPGaveFFTClean(i,:)/std(W21PPGaveFFTClean(i,:))); 

    PPGaveFFTFIN(i,:) = W1PPGaveFFTClean(i,:)+ W2PPGaveFFTClean(i,:) ; % ensambling W1 & W2
end