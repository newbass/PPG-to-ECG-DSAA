function [BPMest, range_Idx] = history_fitting(Idnb, i, BPMest, range_Idx, PPGaveFFTFIN, FreqRange, FreqRangePPG)
    histint = 25;
    
    % History tracking
    if i>30
        
        histint = max(abs(diff(BPMest(1:i-1))))+5;
        
    end

    % HR estimation
    
    if isempty(range_Idx)
        
        [extra, id_x]= max(PPGaveFFTFIN(i,:));
        
        BPMest(i) = FreqRangePPG(id_x(1))*60; 
        
        range_Idx = id_x(1)-round(histint/((FreqRange(2)-FreqRange(1))*60)):id_x(1)+round(histint/((FreqRange(2)-FreqRange(1))*60));
        
    else
        
        [extra, id_x]= max(PPGaveFFTFIN(i,range_Idx));
        
        BPMest(i) = FreqRangePPG(range_Idx(id_x(1)))*60; 
        
        range_Idx = range_Idx(id_x(1))-round(histint/((FreqRange(2)-FreqRange(1))*60)):range_Idx(id_x(1))+round(histint/((FreqRange(2)-FreqRange(1))*60));
        
    end
    
    range_Idx(range_Idx<1) = []; range_Idx(range_Idx>length(FreqRange)) = [];

    % Mild smoothing with linear regression prediction
    
    if i>5 && abs(BPMest(i)-BPMest(i-1))>5
        
        ddd= polyfit(1:length(BPMest(max(1,i-5):i-1)),BPMest(max(1,i-5):i-1),1);
        
        BPMest(i) = 0.8*BPMest(i)+0.2*polyval(ddd,length(BPMest(max(1,i-5):i-1))+1);
        
    end

    %Correction for delay
    mul=0.1;
    
    BPMest(i) = BPMest(i)+sum(sign(BPMest(max(2,i-6):i)-BPMest(max(1,i-7):i-1))*mul);
end