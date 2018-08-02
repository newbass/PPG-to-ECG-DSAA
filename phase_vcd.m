function [FreqRangePPG] = phase_vcd(FreqRangePPG, FreqRange, PPG_ave_FFT, PPG_ave_FFTpr)
    % start phase vocoder for current and previous frames
    for ii=1:size(FreqRangePPG,2)
        curPhase = angle(PPG_ave_FFT(ii));
        prevPhase = angle(PPG_ave_FFTpr(ii)); vocoder = zeros(1,20);
        for n = 1:20
            vocoder(n) = ((curPhase-prevPhase)+(2*pi*(n-1)))/(2*pi*2);
        end
        difference = vocoder - FreqRange(ii);
        [extra, deltaidx] = min(abs(difference));
        FreqRangePPG(ii) = vocoder(deltaidx);
    end
end