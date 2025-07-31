% Model Training
clearvars -except S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 
clc

subjects = {'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'S10'};
leg_side = ['both']; % use this as default if data is missing

subject_weights     = [59.1,59.7,84.1,87.9,70.4,72.8,72.2,89.9,57.4,77.3]; % [kg]
subject_shoe_size   = [42, 40, 42, 42, 42, 40, 40, 42, 40, 42]; % [EU Size]
subject_shoe_size_m = (subject_shoe_size + 10)./200; % [m]

if ~exist('S1','var')
    load S1_S10_AddedInfo.mat;
end

load butterworth_3rd_FS100Hz_FC15Hz.mat; % butterworth filter, 3rd order, fs 100Hz, fc 15Hz

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

% Normalising parameters between [-0.9, 0.9]
max_thumb   =50*ones(1,6);
max_little  =50*ones(1,6);
max_heel    =50*ones(1,6);
max_accel   =[20 20 20]; % [g, with gravity]
max_gyro    =600; % [deg/s]
max_weight  =100; % [N]
max_speed   =1.5; % [m/s]
max_time    =20; % [s]
max_inc     =11; % [deg]
max_shoe    =43; % [EU Size]
max_shoe_mm =0.26; % [m]


%% generate dataset
window=60; % approximately one cycle (data window length)

step_between_points=window/10;   % sliding window
x_train =[]; % for model training
y_train =[];
x_test  =[]; % for model validation after each epoch
y_test  =[];
x_val   =[]; % for model final test
y_val   =[]; 
steps_to_skip_start =300; % neglect the first 300 points after the first HS
steps_to_skip_end   =5; % neglect the last 5 points

total_subject = 1:10;
val_subjects  = 1;
train_test_subjects = [2, 3, 4, 5, 6, 7, 8, 9, 10]; % first 80% for traning and last 20% for testing

% ↓ loop each subject
for i=total_subject
    % ↓ get all fieldnames under 'ShokacData'
    experiments=fieldnames(results.(subjects{i}).ShokacData);

    % ↓ loop each condition
    for j=1:length(experiments) 
        % ↓ get a condition fieldname
        experiment_name=experiments{j};
        % ↓ get speed value [m/s]
        speed=str2num(strrep(experiment_name(strfind(experiment_name,'p')+1:strfind(experiment_name,'p')+3),'x','.'));
        % ↓ get incline angle [deg]
        inc=sscanf(strrep(experiment_name(strfind(experiment_name,'n')+1:end),'_','-'), '%f');

        if isempty(inc)
            break
        end
        % temporarily negect toe clearance collected with non-zero treadmill inclination
        % if inc ~= 0 
        %     continue;
        % end

        % ↓ get all data of this condition under 'ShokacData'
        experiment_data=results.(subjects{i}).ShokacData.(experiment_name);
        % ↓ get all sub-fieldnames of this condition under 'ShokacData'
        sensor_names=fieldnames(experiment_data);
        % ↓ index w.r.t the entry of the last data window
        last_entry=length(experiment_data.time) - steps_to_skip_end - window;
        % ↓ index w.r.t the entry of the first data window (time w.r.t to the first HS + skip steps)
        start_entry = find(experiment_data.time - experiment_data.time(1) ==...
                           experiment_data.shokac_starting_time) + steps_to_skip_start;
        
        
        % ↓ loop each data window with a moving step of 'step_between_points'
        for k = start_entry:step_between_points:last_entry

            % ↓ get synchronized vicon data with the shoe
            eq_vicon_time = experiment_data.time(k:k+window) -...
                            experiment_data.time(1) - experiment_data.time_diff;
            
            % ↓ in case the shoe collection time longer than the vicon
            if max(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).data.Timesteps) < max(eq_vicon_time)
                break;
            end

            % ↓ interpolation for data alignment
            tc = interp1(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).data.Timesteps,...
                         results.(subjects{i}).MarkerData.([experiment_name '_Processed']).toe_clearance_l,...
                         eq_vicon_time);

            % ↓ output/goal (toe distance against the ground) [mm]
            tc_f = filtfilt(SOS, G, tc);
            y_set = {tc_f'};


            % ↓ 3D orientation angle estimated using Marker data
            eulerAgl = interp1(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).data.Timesteps,...
                               results.(subjects{i}).MarkerData.([experiment_name '_Processed']).data.EulerAngles_l,...
                               eq_vicon_time);


            % ↓ create a window of input data (only left)
            x_set = {[  (experiment_data.time(k:k+window,1)-experiment_data.time(k,1))';...

                        % (experiment_data.thumb_l(k:k+window,  2:7)  ./max_thumb)';...
                        % (experiment_data.little_l(k:k+window, 2:7)  ./max_little)';...
                        % (experiment_data.heel_l(k:k+window,   2:7)  ./max_heel)';...

                        % logic grf values result in better accuracy
                        % (capture gait state)
                        (experiment_data.thumb_l(k:k+window,  2:7) > 0.5)';...
                        (experiment_data.little_l(k:k+window, 2:7) > 0.5)';...
                        (experiment_data.heel_l(k:k+window,   2:7) > 0.5)';...
                        
                        % Grf normalised values (may be able to capture -ve
                        % toe clearance due to shoe deformation?)
                        (experiment_data.thumb_l(k:k+window,  2:7))';...
                        (experiment_data.little_l(k:k+window, 2:7))';...
                        (experiment_data.heel_l(k:k+window,   2:7))';...
                        
                        % (experiment_data.accel_l(k:k+window, 1:3)  ./max_accel)';...
                        % (experiment_data.accel_l(k:k+window, 1:3))';...
                        (experiment_data.accel_l(k:k+window, 1:3))'.*9.8;...

                        % (experiment_data.gyro_l(k:k+window, 1)     ./max_gyro)';...
                        (experiment_data.gyro_l(k:k+window, 1))';...

                        % (inc*ones(length(k:k+window),1)/max_inc)';...
                        (inc*ones(length(k:k+window),1))';...

                        eulerAgl';...

                        % (subject_weights(i)*ones(length(k:k+window),1) ./max_weight)';...
                        % (subject_shoe_size(i)*ones(length(k:k+window),1) ./max_shoe)';...
                        (subject_shoe_size_m(i)*ones(length(k:k+window),1))',                      
                        % (experiment_data.phase_l(k:k+window, 1))'
                    ]};


            % fill the dataset 
            if ismember(i, train_test_subjects) % train dataset
                if k < start_entry+(last_entry-start_entry)*0.8
                    x_train = [x_train; x_set];
                    y_train = [y_train; y_set];
                else % test dataset used after each epoch
                    x_test = [x_test; x_set];
                    y_test = [y_test; y_set];
                end
            elseif ismember(i, val_subjects) % validation dataset used finally to test the model
                x_val = [x_val; x_set];
                y_val = [y_val; y_set];
            end

        end

        % display the current progress
        disp([char(subjects{i}) ': ' char(experiments{j})]);

    end
end


%% ↓ remove NaN elements, which may affects the model training
emptyIdx = false(numel(y_train), 1);
for i = 1:numel(y_train)
    yi = y_train{i};
    xi = x_train{i};

    if isempty(yi) || isempty(xi)
        emptyIdx(i) = true;
        continue
    end

    validIdx = ~any(isnan(yi), 1);
    if all(~validIdx)
        emptyIdx(i) = true;
        continue
    end

    y_train{i} = yi(:, validIdx);
    x_train{i} = xi(:, validIdx);
end
y_train(emptyIdx) = [];
x_train(emptyIdx) = [];

% ↓ remove NaN elements
emptyIdx = false(numel(y_test),1);
for i = 1:numel(y_test)
    yi = y_test{i};          
    xi = x_test{i}; 

    if isempty(yi) || isempty(xi)
        emptyIdx(i) = true;
        continue
    end

    validIdx = ~any(isnan(yi), 1);
    if all(~validIdx)
        emptyIdx(i) = true;
        continue
    end

    y_test{i} = yi(:, validIdx);
    x_test{i} = xi(:, validIdx);
end
y_test(emptyIdx) = [];
x_test(emptyIdx) = [];

% ↓ remove NaN elements
emptyIdx = false(numel(y_val),1);
for i = 1:numel(y_val)
    yi = y_val{i};          
    xi = x_val{i}; 

    if isempty(yi) || isempty(xi)
        emptyIdx(i) = true;
        continue
    end

    validIdx = ~any(isnan(yi), 1);
    if all(~validIdx)
        emptyIdx(i) = true;
        continue
    end
    y_val{i} = yi(:, validIdx);
    x_val{i} = xi(:, validIdx);
end
y_val(emptyIdx) = [];
x_val(emptyIdx) = [];

% ↓ scale down the output
scale_y = 0.1;
scale_y = 0.02;
y_train = cellfun(@(v) v * scale_y, y_train, 'UniformOutput', false);
y_test  = cellfun(@(v) v * scale_y, y_test,  'UniformOutput', false);
y_val   = cellfun(@(v) v * scale_y, y_val,   'UniformOutput', false);

%% LSTMs

% ↓ feature dimension of each frame of input data 
numFeatures = size(x_train{1}, 1);
% ↓ two hidden layers
numHiddenUnits1 = numFeatures*6;
numHiddenUnits2 = numFeatures*3;
% ↓ feature dimension of each frame of output data 
numResponses = size(y_train{1}, 1); 
% ↓ para initialization
dropoutLayer1=0.2;
NEpoch=80;
LearningRate=0.002;
% shuffle = 'once';
shuffle = 'every-epoch'; % sequence-to-sequence training and estimation
% layer1='last';
layer1='sequence';
mini_batch_size=128;

if dropoutLayer1>0
    % layers = [  sequenceInputLayer(numFeatures) % input layer (feature dimensions)
    %             lstmLayer(numHiddenUnits1, 'OutputMode', 'sequence')  % 1st layer
    %             fullyConnectedLayer(64)
    %             batchNormalizationLayer
    %             dropoutLayer(dropoutLayer1)
    % 
    %             lstmLayer(numHiddenUnits2, 'OutputMode', 'sequence') % 2nd layer
    %             fullyConnectedLayer(numResponses) % output layer, feature dimension
    %             regressionLayer
    %          ];
    
        layers = [  sequenceInputLayer(numFeatures) % input layer
                bilstmLayer(numFeatures*5, 'OutputMode', 'sequence')  % 1st layer
                dropoutLayer(dropoutLayer1)
                 
                lstmLayer(numFeatures*4, 'OutputMode', 'sequence')  % 2nd layer
                batchNormalizationLayer
                dropoutLayer(dropoutLayer1)

%                 lstmLayer(numFeatures*2, 'OutputMode', 'sequence') % 3rd layer
                lstmLayer(numFeatures*2, 'OutputMode', 'sequence') % 3rd layer
                fullyConnectedLayer(64)
                batchNormalizationLayer
                dropoutLayer(dropoutLayer1)
% 
                fullyConnectedLayer(numResponses) % output layer
                regressionLayer
             ];
else % a simple one
    layers = [  sequenceInputLayer(numFeatures)
                lstmLayer(numHiddenUnits1,OutputMode="sequence")
                fullyConnectedLayer(numResponses)
                regressionLayer];
end

options = trainingOptions('adam',...
                          'ValidationData',{x_test, y_test},...
                          'ValidationPatience', 10,... % times for end without changes
                          'MiniBatchSize', mini_batch_size,... % sample batch size for each training
                          'MaxEpochs', NEpoch, ... % max epoch number
                          'InitialLearnRate', LearningRate, ... % initial learn rate
                          'LearnRateSchedule', 'piecewise', ... % decrease method
                          'LearnRateDropPeriod', 5, ... % decrease rate once each 5 epoch
                          'LearnRateDropFactor', 0.5, ... % decrease 50%
                          'L2Regularization', 1e-4, ...
                          'Shuffle', shuffle, ...
                          'Verbose', false, ...
                          'Plots', 'training-progress'); % pop training progress figure


[net, info] = trainNetwork(x_train, y_train, layers, options);



%% calculate RMSE, R2, Bias, 95%CI

% scale output back
y_true_train = cellfun(@(v) v ./ scale_y, y_train, 'UniformOutput', false);
y_pred_train = cellfun(@(v) v ./ scale_y, predict(net, x_train), 'UniformOutput', false);
y_true_val = cellfun(@(v) v ./ scale_y, y_val, 'UniformOutput', false);
y_pred_val = cellfun(@(v) v ./ scale_y, predict(net, x_val), 'UniformOutput', false);

% [TRAIN DATASET] calculate RMSE
train_mse  = cellfun(@(y_true, y_pred) mean((y_true - y_pred).^2, 'all'), y_true_train, y_pred_train);
train_rmse = sqrt(mean(train_mse));
fprintf('Train RMSE = %.4f mm\n', train_rmse);

% [TRAIN DATASET] calculate R2, Bias, and 95%CI
all_true = cell2mat(y_true_train(:));
all_pred = cell2mat(y_pred_train(:));

y_true_flat = all_true(:);
y_pred_flat = all_pred(:);

SS_res = sum((y_true_flat - y_pred_flat).^2);
SS_tot = sum((y_true_flat - mean(y_true_flat)).^2);

R2_train = 1 - SS_res / SS_tot;

diffs = all_pred - all_true;
diffs_flat = diffs(:); 

bias_train = mean(diffs_flat);
std_dev = std(diffs_flat);
LoA_upper_train = bias_train + 1.96 * std_dev;
LoA_lower_train = bias_train - 1.96 * std_dev;
fprintf('Train R2 = %.4f mm\n', R2_train);
fprintf('Train Bias = %.2f mm\n', bias_train);
fprintf('Train 95CI = [%.2f mm, %.2f mm]\n', LoA_lower_train, LoA_upper_train);

%[VALIDATION DATASET] calculate RMSE
val_mse  = cellfun(@(y_true, y_pred) mean((y_true - y_pred).^2, 'all'), y_true_val, y_pred_val);
val_rmse = sqrt(mean(val_mse));
fprintf('Validation RMSE = %.4f mm\n', val_rmse);

% [VALIDATION DATASET] calculate R2, Bias, and 95%CI
all_true = cell2mat(y_true_val(:));
all_pred = cell2mat(y_pred_val(:));

y_true_flat = all_true(:);
y_pred_flat = all_pred(:);

SS_res = sum((y_true_flat - y_pred_flat).^2);              % Residual sum of squares
SS_tot = sum((y_true_flat - mean(y_true_flat)).^2);        % Total sum of squares

R2_val = 1 - SS_res / SS_tot;

diffs = all_pred - all_true;
diffs_flat = diffs(:); 

bias_val = mean(diffs_flat);
std_dev = std(diffs_flat);
LoA_upper_val = bias_val + 1.96 * std_dev;
LoA_lower_val = bias_val - 1.96 * std_dev;
fprintf('Validation R2 = %.4f mm\n', R2_val);
fprintf('Validation Bias = %.2f mm\n', bias_val);
fprintf('Validation 95CI = [%.2f mm, %.2f mm]\n', LoA_lower_val, LoA_upper_val);


%% plot figures
fontsize    = 15;
linewidth   = 1.5;
blue_colour     = [0.2200 0.6447 0.8410];
orange_colour   = [0.8500 0.3250 0.0980];
yellow_colour   = [0.9290 0.6940 0.1250];


%% train dataset
yy_true_train=[];
for i = 1:length(y_true_train)
    data = y_true_train{i};
    if length(data) < window-step_between_points
        continue;
    end
    yy_true_train = [yy_true_train, data(length(data)-step_between_points+1:end)];
end

yy_pred_train=[];
for i = 1:length(y_pred_train)
    data = y_pred_train{i};
    if length(data) < window-step_between_points
        continue;
    end
    yy_pred_train = [yy_pred_train, data(length(data)-step_between_points+1:end)];
end


minVal = min( [yy_true_train, yy_pred_train] );
maxVal = max( [yy_true_train, yy_pred_train] );

% for extended scope
minVal = minVal-abs(maxVal- minVal)*0.1;
maxVal = maxVal+abs(maxVal- minVal)*0.1;

train_corr_fig=figure; hold on; 
xlabel('yy_true_train (mm)', 'Interpreter', 'none'); 
ylabel('yy_pred_train (mm)', 'Interpreter', 'none');
plot(minVal:maxVal, minVal:maxVal, '--', 'color', 'black', 'LineWidth', 3); % reference line
plot(yy_true_train, yy_pred_train, 'o', 'MarkerEdgeColor', blue_colour);

train_comp_fig=figure; hold on; 
xlabel('points', 'Interpreter', 'none'); 
ylabel('yy_true_train & yy_pred_train (mm)', 'Interpreter', 'none');
plot(1:length(yy_true_train), yy_true_train, 'Color', blue_colour, 'LineWidth', linewidth);
plot(1:length(yy_pred_train), yy_pred_train, '-.', 'Color', orange_colour, 'LineWidth', linewidth);
legend('yy_true_train', 'yy_pred_train', 'Interpreter', 'none');


%% validation dataset
yy_true_val=[];
for i = 1:length(y_true_val)
    data = y_true_val{i};
    if length(data) < window-step_between_points
        continue;
    end
    yy_true_val = [yy_true_val, data(length(data)-step_between_points+1:end)];
end

yy_pred_val=[];
for i = 1:length(y_pred_val)
    data = y_pred_val{i};
    if length(data) < window-step_between_points
        continue;
    end
    yy_pred_val = [yy_pred_val, data(length(data)-step_between_points+1:end)];
end

minVal = min( [yy_true_val, yy_pred_val] );
maxVal = max( [yy_true_val, yy_pred_val] );

% for extended scope
minVal = minVal-abs(maxVal- minVal)*0.1;
maxVal = maxVal+abs(maxVal- minVal)*0.1;

val_corr_fig=figure; hold on; 
xlabel('yy_true_val (mm)', 'Interpreter', 'none'); 
ylabel('yy_pred_val (mm)', 'Interpreter', 'none');
plot(minVal:maxVal, minVal:maxVal, '--', 'color', 'black', 'LineWidth', 3); % reference line
plot(yy_true_val, yy_pred_val, 'o', 'MarkerEdgeColor', blue_colour);

val_comp_fig=figure; hold on; 
xlabel('points', 'Interpreter', 'none'); 
ylabel('yy_true_val & yy_pred_val (mm)', 'Interpreter', 'none');
plot(1:length(yy_true_val), yy_true_val, 'Color', blue_colour, 'LineWidth', linewidth);
plot(1:length(yy_pred_val), yy_pred_val, '-.', 'Color', orange_colour, 'LineWidth', linewidth);
legend('yy_true_val', 'yy_pred_val', 'Interpreter', 'none');

%% Save figs and results
% Define results path
local_data_dir=fullfile(pwd,'..','..','..','Moonshot','LocalData','Experiment','FESTrajectories_v02','TrainingModelsResults');
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

% Save figures and useful script information 
filename='TrainCorrFig.png';
saveas(train_corr_fig,[local_data_dir filesep filename])
filename='TrainCorrFig.fig';
saveas(train_corr_fig,[local_data_dir filesep filename])

filename='TrainCompFig.png';
saveas(train_comp_fig,[local_data_dir filesep filename])
filename='TrainCompFig.fig';
saveas(train_comp_fig,[local_data_dir filesep filename])

filename='ValCorrFig.png';
saveas(val_corr_fig,[local_data_dir filesep filename])
filename='ValCorrFig.fig';
saveas(val_corr_fig,[local_data_dir filesep filename])

filename='ValCompFig.png';
saveas(val_comp_fig,[local_data_dir filesep filename])
filename='ValCompFig.fig';
saveas(val_comp_fig,[local_data_dir filesep filename])

ModelTuningParams = struct('dropoutLayer1',dropoutLayer1,'LearningRate',LearningRate,'mini_batch_size',mini_batch_size,'NEpoch',NEpoch,'numFeatures',numFeatures,'numHiddenUnits1',numHiddenUnits1,'numHiddenUnits2',numHiddenUnits2,'numResponses',numResponses,'options',options,'scale_y',scale_y,'window',window,'x_set',x_set);
ResultsParams = struct('bias_train',bias_train,'bias_val',bias_val,'LoA_lower_train',LoA_lower_train,'LoA_upper_train',LoA_upper_train,'LoA_lower_val',LoA_lower_val,'LoA_upper_val',LoA_upper_val,'R2_train',R2_train,'R2_val',R2_val,'train_rmse',train_rmse,'train_test_subjects',train_test_subjects,'val_rmse',val_rmse,'val_subjects',val_subjects);

save([local_data_dir filesep 'Trial_details.mat'],'ModelTuningParams','ResultsParams')