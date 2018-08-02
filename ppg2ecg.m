% Authors: Andriy Temko, INFANT Research Centre
% http://www.infantcentre.ie/
% http://eleceng.ucc.ie/~andreyt/
function ppg2ecg(path_to_train, path_to_test, path_to_save)
    
    % Norm for path
    path_to_train = char(path_to_train);
    if(path_to_train(end) ~= '/')
        path_to_train(end + 1) = '/';
    end 
    path_to_test = char(path_to_test);
    if(path_to_test(end) ~= '/')
        path_to_test(end + 1) = '/';
    end
    path_to_save = char(path_to_save);
    if(path_to_save(end) ~= '/')
        path_to_save(end + 1) = '/';
    end
    
    % Getting all files
    train_files = dir(char(string(path_to_train) + "*.mat"));
    test_files = dir(char(string(path_to_test) + "*.mat"));

    train_files_length = 0;
    ALData = {};
    ECData = {};

    % Storing all file names
    for i = 1:length(train_files)
        if(~length(regexp(train_files(i).name, "BPMtrace")))
            ALData{end + 1} = train_files(i).name;
            train_files_length = train_files_length + 1;
        else
            ECData{end + 1} = train_files(i).name;
        end
    end
    for i = 1:length(test_files)
        if(~length(regexp(test_files(i).name, "BPMtrace")))
            ALData{end + 1} = test_files(i).name;
        end
    end
    
    % Overall parameters
    [srate, FFTres, WFlength, total_count] = deal(125, 1024, 10, size(ALData,2));

    % Filter parameters
    [Low_freq, High_freq] = deal(0.5, 3);
    [b_coeff,a_coeff] = butter(4, [0.4 4]/(125/2),'bandpass');

    % Metrics
    fullBPM=[];fullBPM0=[]; % for computing the correlation
    [my_Error, my_Error_Std, my_Rel_Error] = deal(zeros(1,total_count), zeros(1,total_count), zeros(1,total_count));

    % framework rule, 8s window 2s shift, no look into future
    window   = 8 * srate;  % window length is 8 seconds
    step     = 2 * srate;  % step size is 2 seconds

    % loop for each file
    for Idnb = 1 : total_count
        % load the data
        if Idnb > train_files_length
            load(char(string(path_to_test) + string(ALData{Idnb})));
            [ch1, ch2, ch3, ch4, ch5] = deal(1, 2, 3, 4, 5);
        else
            load(char(string(path_to_train) + string(ALData{Idnb})));
            [ch1, ch2, ch3, ch4, ch5] = deal(2, 3, 4, 5, 6);
        end
        
        windowNb = floor((length(sig)-window)/step) + 1;  % total number of windows(estimates)

        % initialization of variables
        BPMest = zeros(1,windowNb); % estimated BPM
        range_Idx = []; % range of search for the next estimates
        clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean W11_PPG_ave_FFT_Clean PPG_ave_FFT_FIN W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN;
        PPGaveFFTFIN = [];
        [W1FFT, W11FFT, W2FFT, W21FFT] = deal([]);
        
        % Traversal through window
        for i =  [1 :  windowNb]
            cur_Segment = (i-1)*step+1 : (i-1)*step+window;
            cur_Data = sig(:,cur_Segment);

            PPG1 = cur_Data(ch1,:); PPG2 = cur_Data(ch2,:);
            ACCX = cur_Data(ch3,:); ACCY = cur_Data(ch4,:); ACCZ = cur_Data(ch5,:);

            % Reducing given data to other standargs
            [PPG1, PPG2, ACCX, ACCY, ACCZ, PPG_ave, srate, FreqRange, PPG_ave_FFT, ACC_X_FFT, ACC_Y_FFT, ACC_Z_FFT] = data_reconfig(PPG1, PPG2, ACCX, ACCY, ACCZ, b_coeff, a_coeff, FFTres, High_freq, Low_freq);

            % phase vocoder to refine spectral estimations
            Freq_Range_PPG = FreqRange;
            if i>1
                [Freq_Range_PPG] = phase_vcd(Freq_Range_PPG, FreqRange, PPG_ave_FFT, PPG_ave_FFTpr);
            end

            % smooth phase vocoder frequency estimates
            Freq_Range_PPG = moving_average(Freq_Range_PPG,3);

            % save previous spectrum for the next phase vocoder call
            PPG_ave_FFTpr = PPG_ave_FFT;

            % Wiener filtering
            [W1FFT, W11FFT, W2FFT, W21FFT, PPGaveFFTFIN] = wiener_filter(i, WFlength, PPG_ave_FFT, ACC_X_FFT, ACC_Y_FFT, ACC_Z_FFT, W1FFT, W11FFT, W2FFT, W21FFT, PPGaveFFTFIN);

            % Linear regression and polygraphic substitution
            [BPMest, range_Idx] = history_fitting(Idnb, i, BPMest, range_Idx, PPGaveFFTFIN, FreqRange, Freq_Range_PPG);
        end

        % Code to store data
        if Idnb > train_files_length
            pred = BPMest.';
            pred = pred(1:125);
            save([path_to_save 'output_team_29_' ALData{Idnb}], 'pred');
        else
            load(char(string(path_to_train) + string(ECData{Idnb})));           % load groundtruth
            [my_Error(Idnb), my_Rel_Error(Idnb), my_Error_Std(Idnb)] = deal(mean(abs(BPM0 - BPMest(1:1:end)')), mean(abs(BPM0 - BPMest(1:1:end)')./BPM0), std(abs(BPM0 - BPMest(1:1:end)')));
        end

    end
    
    % Error display
    fprintf('Err12=%2.2f(%2.2f), Err11=%2.2f(%2.2f), ErrAll=%2.2f(%2.2f)\n', mean(my_Error(1:12)),mean(my_Error_Std(1:12)),mean(my_Error(13:end)),mean(my_Error_Std(13:end)),mean(my_Error),mean(my_Error_Std));
    for s=1:total_count
        fprintf(' %2.2f', my_Error(s));
    end
    fprintf('\n');
    [max_val, max_index] = max(my_Error(1:12));
    [min_val, min_index] = min(my_Error(1:12));
    fprintf("Maximum error occurs at %d and minimum error occurs at %d\n", max_index, min_index);
end