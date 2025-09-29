total_subject = 1:10;

%% leave-one-subject-out
for testSubject = total_subject
    subjects = {'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'S10'};
    subject_shoe_size = [42, 40, 42, 42, 42, 40, 40, 42, 40, 42]; % [EU Size]

    load butterworth_3rd_FS100Hz_FC15Hz.mat; % butterworth filter, 3rd order, fs 100Hz, fc 15Hz
    load('S1_S10_AddedInfo_05-Aug.mat');

    results.S1  =S1;
    results.S2  =S2;
    results.S3  =S3;
    results.S4  =S4;
    results.S5  =S5;
    results.S6  =S6;
    results.S7  =S7;
    results.S8  =S8;
    results.S9  =S9;
    results.S10 =S10;

    % datasets for model training
    All_Output.tc    = [];
    All_Input.shoeTime = [];
    All_Input.thumb  = [];
    All_Input.little = [];
    All_Input.heel   = [];
    All_Input.accel  = [];
    All_Input.gyro   = [];
    All_Input.inc      = [];
    All_Input.shoeSize = [];
    All_Input.legSide  = [];
    All_Input.EulerAngles_footIMU = [];

    % datasets for model test
    All_OutputTes.tc    = [];
    All_InputTes.shoeTime = [];
    All_InputTes.thumb  = [];
    All_InputTes.little = [];
    All_InputTes.heel   = [];
    All_InputTes.accel  = [];
    All_InputTes.gyro   = [];
    All_InputTes.inc      = [];
    All_InputTes.shoeSize = [];
    All_InputTes.legSide  = [];
    All_InputTes.EulerAngles_footIMU = [];

    legSide = ["left", "right"];
    edgeCase = 'exclu'; % inclu, exclu, only (only process the edge case data)

    % combine data across subjects and conditions
    for side = legSide
        for i = total_subject
            steps_to_skip_start =300; % neglect the first 300 points after the first HS
            steps_to_skip_end   =5; % neglect the last 5 points

            experiments =fieldnames(results.(subjects{i}).ShokacData);

            % ↓ loop each condition
            for j=1:length(experiments)
                experiment_name=experiments{j};
                inc=sscanf(strrep(experiment_name(strfind(experiment_name,'n')+1:end),'_','-'), '%f');

                % edge case data
                if strcmp(edgeCase, 'exclu')
                    if contains(experiment_name, 'edge')
                        continue;
                    end
                elseif strcmp(edgeCase, 'inclu')
                elseif strcmp(edgeCase, 'only')
                    if ~contains(experiment_name, 'edge')
                        continue;
                    end
                else
                    warning('error in edgeCase');
                    pause;
                end

                % get all data of this condition under 'ShokacData'
                experiment_data=results.(subjects{i}).ShokacData.(experiment_name);
                endPoint    = length(experiment_data.time) - steps_to_skip_end;
                startPoint  = find(experiment_data.time - experiment_data.time(1) ==...
                    experiment_data.shokac_starting_time) + steps_to_skip_start;

                index = startPoint:endPoint;
                eq_vicon_time = experiment_data.time(index) - experiment_data.time(1) - experiment_data.time_diff;

                % if max Vicon time is less than selected shoe time, then reduce
                % the shoe time index
                maxViconTime = max(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).Timesteps);
                if maxViconTime < max(eq_vicon_time)
                    indd = find((eq_vicon_time - maxViconTime) < 0);
                    indexDiff = length(eq_vicon_time) - indd(end);
                    index = index(1: length(eq_vicon_time) - indexDiff);
                end

                % remove the lines from index if there are corrupted elements
                if strcmp(side, 'left')
                    sub_data = experiment_data.masked_data(index, 1:5);
                elseif strcmp(side, 'right')
                    sub_data = experiment_data.masked_data(index, 6:10);
                else
                    warning('error in legSide');
                    pause;
                end
                rows_to_keep = all(sub_data == 1, 2);
                index = index(rows_to_keep);

                % if length(rows_to_keep) ~= sum(rows_to_keep)
                %     disp('There are corrupted lines');
                % end

                % get the time using the latest index
                eq_vicon_time = experiment_data.time(index) - experiment_data.time(1) - experiment_data.time_diff;

                % output: toe clearance
                if strcmp(side, 'left')
                    toe_clearance = results.(subjects{i}).MarkerData.([experiment_name '_Processed']).toe_clearance_l;
                    if i == testSubject
                        All_InputTes.legSide = [All_InputTes.legSide; zeros(length(index), 1)]; % left index 0
                    else
                        All_Input.legSide = [All_Input.legSide; zeros(length(index), 1)]; % left index 0
                    end
                elseif strcmp(side, 'right')
                    toe_clearance = results.(subjects{i}).MarkerData.([experiment_name '_Processed']).toe_clearance_r;
                    if i == testSubject
                        All_InputTes.legSide = [All_InputTes.legSide; ones(length(index), 1)]; % right index 1
                    else
                        All_Input.legSide = [All_Input.legSide; ones(length(index), 1)]; % right index 1
                    end
                else
                    warning('error in legSide');
                    pause;
                end
                tc = interp1(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).Timesteps,...
                    toe_clearance, eq_vicon_time);
                if i == testSubject
                    All_OutputTes.tc = [All_OutputTes.tc; filtfilt(SOS, G, tc)];
                else
                    All_Output.tc = [All_Output.tc; filtfilt(SOS, G, tc)];
                end
                % input: EulerAngles_footIMU
                if strcmp(side, 'left')
                    footIMU = results.(subjects{i}).MarkerData.([experiment_name '_Processed']).EulerAngles_footIMU_l;
                elseif strcmp(side, 'right')
                    footIMU = results.(subjects{i}).MarkerData.([experiment_name '_Processed']).EulerAngles_footIMU_r;
                else
                    warning('error in legSide');
                    pause;
                end
                EulerAngles_footIMU = ...
                    interp1(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).Timesteps,...
                    footIMU', eq_vicon_time);

                if i == testSubject
                    All_InputTes.EulerAngles_footIMU = [All_InputTes.EulerAngles_footIMU; EulerAngles_footIMU];
                else
                    All_Input.EulerAngles_footIMU = [All_Input.EulerAngles_footIMU; EulerAngles_footIMU];
                end
                % input: time
                timeTemp = experiment_data.time(index);
                time = timeTemp-timeTemp(1);
                % time = experiment_data.time(2:end)-experiment_data.time(1:end-1);
                % time = [time;time(end)];

                if i == testSubject
                    All_InputTes.shoeTime = [All_InputTes.shoeTime; time];
                else
                    All_Input.shoeTime = [All_Input.shoeTime; time];
                end

                if strcmp(side, 'left')
                    if i == testSubject
                        % input: force
                        All_InputTes.thumb  = [All_InputTes.thumb; experiment_data.thumb_l(index,  2:7)];
                        All_InputTes.little = [All_InputTes.little; experiment_data.little_l(index,  2:7)];
                        All_InputTes.heel   = [All_InputTes.heel; experiment_data.heel_l(index,  2:7)];

                        % input: acc
                        All_InputTes.accel  = [All_InputTes.accel; experiment_data.accel_l(index, 1:3)];

                        % input: gyro
                        All_InputTes.gyro   = [All_InputTes.gyro; experiment_data.gyro_l(index, 1)];

                        % input: inc
                        All_InputTes.inc      = [All_InputTes.inc; inc*ones( length(experiment_data.gyro_l(index, 1)), 1)];

                        % input: shoe size
                        All_InputTes.shoeSize = [All_InputTes.shoeSize; subject_shoe_size(i)*ones( length(experiment_data.gyro_l(index, 1)), 1)];

                    else
                        % input: force
                        All_Input.thumb  = [All_Input.thumb; experiment_data.thumb_l(index,  2:7)];
                        All_Input.little = [All_Input.little; experiment_data.little_l(index,  2:7)];
                        All_Input.heel   = [All_Input.heel; experiment_data.heel_l(index,  2:7)];

                        % input: acc
                        All_Input.accel  = [All_Input.accel; experiment_data.accel_l(index, 1:3)];

                        % input: gyro
                        All_Input.gyro   = [All_Input.gyro; experiment_data.gyro_l(index, 1)];

                        % input: inc
                        All_Input.inc      = [All_Input.inc; inc*ones( length(experiment_data.gyro_l(index, 1)), 1)];

                        % input: shoe size
                        All_Input.shoeSize = [All_Input.shoeSize; subject_shoe_size(i)*ones( length(experiment_data.gyro_l(index, 1)), 1)];

                    end
                elseif strcmp(side, 'right')
                    if i == testSubject
                        % input: force
                        All_InputTes.thumb  = [All_InputTes.thumb; experiment_data.thumb_r(index,  2:7)];
                        All_InputTes.little = [All_InputTes.little; experiment_data.little_r(index,  2:7)];
                        All_InputTes.heel   = [All_InputTes.heel; experiment_data.heel_r(index,  2:7)];

                        % input: acc
                        All_InputTes.accel  = [All_InputTes.accel; experiment_data.accel_r(index, 1:3)];

                        % input: gyro
                        All_InputTes.gyro   = [All_InputTes.gyro; experiment_data.gyro_r(index, 1)];

                        % input: inc
                        All_InputTes.inc      = [All_InputTes.inc; inc*ones( length(experiment_data.gyro_r(index, 1)), 1)];

                        % input: shoe size
                        All_InputTes.shoeSize = [All_InputTes.shoeSize; subject_shoe_size(i)*ones( length(experiment_data.gyro_r(index, 1)), 1)];
                    else
                        % input: force
                        All_Input.thumb  = [All_Input.thumb; experiment_data.thumb_r(index,  2:7)];
                        All_Input.little = [All_Input.little; experiment_data.little_r(index,  2:7)];
                        All_Input.heel   = [All_Input.heel; experiment_data.heel_r(index,  2:7)];

                        % input: acc
                        All_Input.accel  = [All_Input.accel; experiment_data.accel_r(index, 1:3)];

                        % input: gyro
                        All_Input.gyro   = [All_Input.gyro; experiment_data.gyro_r(index, 1)];

                        % input: inc
                        All_Input.inc      = [All_Input.inc; inc*ones( length(experiment_data.gyro_r(index, 1)), 1)];

                        % input: shoe size
                        All_Input.shoeSize = [All_Input.shoeSize; subject_shoe_size(i)*ones( length(experiment_data.gyro_r(index, 1)), 1)];
                    end
                else
                    warning('error in legSide');
                    pause;
                end

                % display progress
                disp(['No.(test): ' num2str(i) '(' num2str(testSubject) '), ' char(side) ' ' char(subjects{i}) ': ' char(experiments{j})]);
            end
        end
    end

    % calculate Mean from training datasets
    allMean_tc       = mean(All_Output.tc);
    allMean_thumb    = mean(All_Input.thumb);
    allMean_little   = mean(All_Input.little);
    allMean_heel     = mean(All_Input.heel);
    allMean_accel    = mean(All_Input.accel);
    allMean_gyro     = mean(All_Input.gyro);
    allMean_inc      = mean(All_Input.inc);
    allMean_shoeSize = mean(All_Input.shoeSize);
    allMean_EulerAngles_footIMU = mean(All_Input.EulerAngles_footIMU);

    % calculate Std from training datasets
    allStd_tc       = std(All_Output.tc);
    allStd_thumb    = std(All_Input.thumb);
    allStd_little   = std(All_Input.little);
    allStd_heel     = std(All_Input.heel);
    allStd_accel    = std(All_Input.accel);
    allStd_gyro     = std(All_Input.gyro);
    allStd_inc      = std(All_Input.inc);
    allStd_shoeSize = std(All_Input.shoeSize);
    allStd_EulerAngles_footIMU = std(All_Input.EulerAngles_footIMU);

    % z-score normalize training datasets
    nAll_Output.tc      = (All_Output.tc      - allMean_tc)     ./allStd_tc;
    nAll_Input.shoeTime = All_Input.shoeTime;
    nAll_Input.thumb    = (All_Input.thumb    - allMean_thumb)  ./allStd_thumb;
    nAll_Input.little   = (All_Input.little   - allMean_little) ./allStd_little;
    nAll_Input.heel     = (All_Input.heel     - allMean_heel)   ./allStd_heel;
    nAll_Input.accel    = (All_Input.accel    - allMean_accel)  ./allStd_accel;
    nAll_Input.gyro     = (All_Input.gyro     - allMean_gyro)   ./allStd_gyro;
    nAll_Input.inc      = (All_Input.inc      - allMean_inc)    ./allStd_inc;
    nAll_Input.shoeSize = (All_Input.shoeSize - allMean_shoeSize) ./allStd_shoeSize;
    nAll_Input.EulerAngles_footIMU = (All_Input.EulerAngles_footIMU - allMean_EulerAngles_footIMU)./allStd_EulerAngles_footIMU;
    nAll_Input.legSide  = All_Input.legSide;

    % z-score normalize test datasets using Mean and Std calculated from training datasets
    nAll_OutputTes.tc      = (All_OutputTes.tc      - allMean_tc)     ./allStd_tc;
    nAll_InputTes.shoeTime = All_InputTes.shoeTime;
    nAll_InputTes.thumb    = (All_InputTes.thumb    - allMean_thumb)  ./allStd_thumb;
    nAll_InputTes.little   = (All_InputTes.little   - allMean_little) ./allStd_little;
    nAll_InputTes.heel     = (All_InputTes.heel     - allMean_heel)   ./allStd_heel;
    nAll_InputTes.accel    = (All_InputTes.accel    - allMean_accel)  ./allStd_accel;
    nAll_InputTes.gyro     = (All_InputTes.gyro     - allMean_gyro)   ./allStd_gyro;
    nAll_InputTes.inc      = (All_InputTes.inc      - allMean_inc)    ./allStd_inc;
    nAll_InputTes.shoeSize = (All_InputTes.shoeSize - allMean_shoeSize) ./allStd_shoeSize;
    nAll_InputTes.EulerAngles_footIMU = (All_InputTes.EulerAngles_footIMU - allMean_EulerAngles_footIMU)./allStd_EulerAngles_footIMU;
    nAll_InputTes.legSide  = All_InputTes.legSide;

    %% calculate average max window duration using training datasets only
    window=60;
    step_between_points=window/10;   % sliding window
    totalLen  = length(nAll_Output.tc);
    maxElement = [];
    for i = 1:step_between_points:totalLen-window
        index = i:i+window-1;
        temp = nAll_Input.shoeTime(index, :);
        if max(temp) - min(temp) > 5
            continue;
        end
        temp = temp - temp(1);

        maxElement = [maxElement; max(temp)];
    end
    averageMaxWindowDuration = mean(maxElement);


    %% Model training
    % keep consistency
    rng(123, 'threefry'); % for randperm
    gpurng(123, 'threefry'); % for model training

    % Generate datasets for model training
    percentTra = 0.8;
    percentVal = 0.2;

    x_total = [];
    y_total = [];
    totalLen  = length(nAll_Output.tc);
    for i = 1:step_between_points:totalLen-window
        % display progress
        index = i:i+window-1;
        disp(['TRAIN - test No.' num2str(testSubject) ', progress: ' num2str(i/(totalLen-window)*100, '%.3f%%')]);

        % output
        y_set = {nAll_Output.tc(index)'};

        % input
        temp = nAll_Input.shoeTime(index, :);
        if max(temp) - min(temp) > 5 % in case of the window crossing two trials
            continue;
        end
        temp = (temp - temp(1)) / averageMaxWindowDuration;

        x_set = {[
            temp'; ...
            nAll_Input.accel(index, :)'; ...
            nAll_Input.EulerAngles_footIMU(index, [1,3])';...
            nAll_Input.gyro(index, :)';...
            nAll_Input.inc(index, :)';...
            nAll_Input.shoeSize(index, :)';...
            nAll_Input.thumb(index, :)';...
            nAll_Input.legSide(index, :)']};

        y_total = [y_total; y_set];
        x_total = [x_total; x_set];
    end
    [x_total, y_total] = dropNaNCells(x_total, y_total);

    % Separate dataset
    lengthCell = length(y_total);
    length_tra = floor(lengthCell * percentTra);
    length_val = floor(lengthCell * percentVal);

    % shuffle
    idx_train = randperm(lengthCell);
    x_total = x_total(idx_train);
    y_total = y_total(idx_train);

    % generate
    x_tra = x_total(1:length_tra);
    y_tra = y_total(1:length_tra);
    x_val = x_total(end-length_val+1:end);
    y_val = y_total(end-length_val+1:end);

    % Generate test datasets
    x_tes = [];
    y_tes = [];
    totalLen  = length(nAll_OutputTes.tc);
    for i = 1:step_between_points:totalLen-window
        % show the percent progress
        index = i:i+window-1;
        disp(['TEST - test No.' num2str(testSubject) ', progress: ' num2str(i/(totalLen-window)*100, '%.3f%%')]);

        % output
        y_set = {nAll_OutputTes.tc(index)'};

        % input
        temp = nAll_InputTes.shoeTime(index, :);
        if max(temp) - min(temp) > 5 % in case of the window crossing two trials
            continue;
        end
        temp = (temp - temp(1)) / averageMaxWindowDuration;

        x_set = {[
            temp'; ...
            nAll_InputTes.accel(index, :)'; ...
            nAll_InputTes.EulerAngles_footIMU(index, [1,3])';...
            nAll_InputTes.gyro(index, :)';...
            nAll_InputTes.inc(index, :)';...
            nAll_InputTes.shoeSize(index, :)';...
            nAll_InputTes.thumb(index, :)';...
            nAll_InputTes.legSide(index, :)']};

        y_tes = [y_tes; y_set];
        x_tes = [x_tes; x_set];
    end
    [x_tes, y_tes] = dropNaNCells(x_tes, y_tes);

    % shuffle test datasets
    idx_test = randperm( length(y_tes) );
    x_tes = x_tes(idx_test);
    y_tes = y_tes(idx_test);

    % LSTMs
    numFeatures = size(x_tra{1}, 1);
    numHiddenUnits1 = numFeatures*8; % [1] after optimization 26-Aug-2025
    numHiddenUnits2 = numFeatures*8; % [2] after optimization 26-Aug-2025
    numResponses = size(y_tra{1}, 1);
    dropoutLayer1=0.05; % [3] after optimization 26-Aug-2025
    NEpoch=80;
    LearningRate=1e-2; % [4] after optimization 26-Aug-2025
    shuffle = 'every-epoch';
    layer1='sequence';
    mini_batch_size=512; % [5] after optimization 26-Aug-2025

    layers = [
        sequenceInputLayer(numFeatures) % input layer
        bilstmLayer(numHiddenUnits1, 'OutputMode', 'sequence')
        batchNormalizationLayer
        geluLayer
        dropoutLayer(dropoutLayer1)

        lstmLayer(numHiddenUnits2, 'OutputMode', 'sequence')
        batchNormalizationLayer
        geluLayer
        dropoutLayer(dropoutLayer1)

        fullyConnectedLayer(numResponses)
        regressionLayer
        ];

    options = trainingOptions('adam',...
        'ValidationData',{x_val, y_val},...
        'ValidationPatience', 20,...
        'MiniBatchSize', mini_batch_size,...
        'MaxEpochs', NEpoch, ...
        'InitialLearnRate', LearningRate, ...
        'LearnRateSchedule', 'piecewise', ...
        'LearnRateDropPeriod', 4, ...
        'LearnRateDropFactor', 0.6, ...
        'GradientThreshold', 1, ...
        'L2Regularization', 1e-3, ... % [6] after optimization 26-Aug-2025
        'Shuffle', shuffle, ...
        'Verbose', false, ...
        'Plots', 'training-progress',...
        'ExecutionEnvironment', 'auto');

    [net, info] = trainNetwork(x_tra, y_tra, layers, options);

    % calculate RMSE, R2, Bias, 95%CI
    fprintf('\nnet train datasets:\n');
    [y_true_tra, y_pred_tra, indicator_tra] = ...
        printIndicator(x_tra, y_tra, allStd_tc, allMean_tc, net);

    fprintf('\nnet test datasets:\n');
    [y_true_tes, y_pred_tes, indicator_tes] = ...
        printIndicator(x_tes, y_tes, allStd_tc, allMean_tc, net);


    %% Save results
    local_data_dir=fullfile(pwd, 'LOSOResults');
    % Observe present folders (named in arithmetic order)
    files = dir(local_data_dir);
    dirFlags= [files.isdir];
    subFolder = files(dirFlags);
    subFolderNames = [str2double({subFolder(3:end).name})];
    % Observe latest folder of results (largest number)
    maxFolderNum=max(subFolderNames);
    % Create the next folder
    nextFolderName = maxFolderNum+1;
    if nextFolderName<10
        nextFolderName=['00' num2str(nextFolderName)];
    elseif nextFolderName<100
        nextFolderName=['0' num2str(nextFolderName)];
    else
        nextFolderName=['0' num2str(nextFolderName)];
    end
    mkdir(fullfile(local_data_dir,nextFolderName))
    local_data_dir=fullfile(local_data_dir,nextFolderName);

    % Save results
    Results = struct('model_net',net,...
        'model_info',info,...
        'indicator_tra',indicator_tra,...
        'indicator_tes',indicator_tes,...
        'length_tra',length_tra,...
        'length_val',length_val,...
        'idx_train',idx_train,...
        'idx_test',idx_test,...
        'allMean_tc',allMean_tc,...
        'allStd_tc',allStd_tc,...
        'averageMaxWindowDuration', averageMaxWindowDuration,...
        'testSubject', testSubject);

    save([local_data_dir filesep 'Results.mat'],'Results');


    % clear variables for next LOSO test
    clearvars -except total_subject
end


%% Functions
function [xOut, yOut] = dropNaNCells(xIn, yIn)
emptyIdx = false(numel(yIn), 1);
for k = 1:numel(yIn)
    yi = yIn{k};
    xi = xIn{k};
    if isempty(yi) || isempty(xi)
        emptyIdx(k) = true; continue
    end
    validIdx = ~any(isnan(yi),1) & ~any(isnan(xi),1);
    if all(~validIdx)
        emptyIdx(k) = true; continue
    end
    yIn{k} = yi(:, validIdx);
    xIn{k} = xi(:, validIdx);
end
yIn(emptyIdx) = [];
xIn(emptyIdx) = [];
xOut = xIn; yOut = yIn;
end

function [y_true, y_pred, indicator] = printIndicator(x_set, y_set, allStd_tc, allMean_tc, net)
% scale output back
y_true = cellfun(@(v) v .* allStd_tc + allMean_tc, y_set, 'UniformOutput', false);
y_pred = cellfun(@(v) v .* allStd_tc + allMean_tc, predict(net, x_set), 'UniformOutput', false);

% calculate RMSE
data_mse  = cellfun(@(y_true, y_pred) mean((y_true - y_pred).^2, 'all'), y_true, y_pred);
data_rmse = sqrt(mean(data_mse));

% calculate R2
all_true = cell2mat(y_true(:));
all_pred = cell2mat(y_pred(:));

y_true_flat = all_true(:);
y_pred_flat = all_pred(:);

SS_res = sum((y_true_flat - y_pred_flat).^2);
SS_tot = sum((y_true_flat - mean(y_true_flat)).^2);
R2 = 1 - SS_res / SS_tot;

% calculate Bias
diffs = all_pred - all_true;
diffs_flat = diffs(:);
bias = mean(diffs_flat);

% calculate 95%CI
std_dev = std(diffs_flat);
LoA_upper = bias + 1.96 * std_dev;
LoA_lower = bias - 1.96 * std_dev;

indicator.rmse = data_rmse;
indicator.R2 = R2;
indicator.bias = bias;
indicator.LoA_lower = LoA_lower;
indicator.LoA_upper = LoA_upper;

% print results
fprintf('RMSE = %.4f mm\n', data_rmse);
fprintf('R2 = %.4f \n', R2);
fprintf('Bias = %.2f mm\n', bias);
fprintf('95CI = [%.2f, %.2f] mm\n', LoA_lower, LoA_upper);
end

function [corr_fig, comp_fig] = plotFigures(y_true, y_pred, window, step_between_points, name)
linewidth   = 1.5;
blue_colour     = [0.2200 0.6447 0.8410];
orange_colour   = [0.8500 0.3250 0.0980];

yy_true=[];
for i = 1:length(y_true)
    data = y_true{i};
    if length(data) < window-step_between_points
        continue;
    end
    yy_true = [yy_true, data(length(data)-step_between_points+1:end)];
end

yy_pred=[];
for i = 1:length(y_pred)
    data = y_pred{i};
    if length(data) < window-step_between_points
        continue;
    end
    yy_pred = [yy_pred, data(length(data)-step_between_points+1:end)];
end

% minVal = min( [yy_true, yy_pred] );
% maxVal = max( [yy_true, yy_pred] );
% minVal = minVal-abs(maxVal- minVal)*0.1;
% maxVal = maxVal+abs(maxVal- minVal)*0.1;
%
% corr_fig=figure; hold on; title(name);
% xlabel('yy_true (mm)', 'Interpreter', 'none');
% ylabel('yy_pred (mm)', 'Interpreter', 'none');
% plot(minVal:maxVal, minVal:maxVal, '--', 'color', 'black', 'LineWidth', 3); % reference line
% plot(yy_true, yy_pred, 'o', 'MarkerEdgeColor', blue_colour);

comp_fig=figure; hold on; title(name);
xlabel('points', 'Interpreter', 'none');
ylabel('yy_true & yy_pred (mm)', 'Interpreter', 'none');
plot(1:length(yy_true), yy_true, 'Color', blue_colour, 'LineWidth', linewidth);
plot(1:length(yy_pred), yy_pred, '-.', 'Color', orange_colour, 'LineWidth', linewidth);
legend('yy_true', 'yy_pred', 'Interpreter', 'none');

end
