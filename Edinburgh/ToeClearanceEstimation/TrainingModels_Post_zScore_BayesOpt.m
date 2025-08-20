% Model Training

percentVal = [0.0, 0.1]; % 0-10% for validation
percentTra = [0.1, 0.8]; % 10-80% for training
percentTes = [0.8, 1.0]; % 80-100% for testing

window=60; % approximately one cycle (data window length) % could be optimised but at a high computational cost

step_between_points=window/10;   % sliding window

if ~exist('x_tra','var')
    
    x_tra =[]; % for model training
    y_tra =[];
    x_tes =[]; % for model test after each epoch
    y_tes =[];
    x_val =[]; % for model final validation
    y_val =[];
    
    totalLen  = length(nAll_Output.tc_l);
    index_val = round(percentVal .* totalLen);
    index_tra = round(percentTra .* totalLen);
    index_tes = round(percentTes .* totalLen);
    
    
    for i = 1:step_between_points:totalLen-window
        
        % show the percent progress
        index = i:i+window-1;
        disp(['progress: ' num2str(i/(totalLen-window)*100, '%.3f%%')]);
        
        % output
        y_set = {nAll_Output.tc_l(index)'};
        
        % input
        x_set = {[
            nAll_Input.shoeTime(index, :)'; ...
            nAll_Input.accel_l(index, :)'; ...
            nAll_Input.EulerAngles_footIMU_l(index, :)';...
            nAll_Input.gyro_l(index, :)';...
            nAll_Input.heel_l(index, :)';...
            nAll_Input.inc(index, :)';...
            nAll_Input.little_l(index, :)';...
            nAll_Input.shoeSize(index, :)';...
            nAll_Input.thumb_l(index, :)']};
        
        % put them into vectors
        if i>=min(index_val) && i<=max(index_val)
            y_val = [y_val; y_set];
            x_val = [x_val; x_set];
        elseif i>=min(index_tra) && i<=max(index_tra)
            y_tra = [y_tra; y_set];
            x_tra = [x_tra; x_set];
        elseif i>=min(index_tes) && i<=max(index_tes)
            y_tes = [y_tes; y_set];
            x_tes = [x_tes; x_set];
        end
        
    end
    
    x_train = x_tra;
    y_train = y_tra;
    x_test  = x_tes;
    y_test  = y_tes;

    % ↓ remove NaN elements, which may affects the model training
    [x_train,y_train] = dropNaNCells(x_train,y_train);
    [x_val,  y_val  ] = dropNaNCells(x_val,  y_val  );
    [x_test, y_test ] = dropNaNCells(x_test, y_test );

end

%% LSTMs

% ↓ feature dimension of each frame of input data
numFeatures = size(x_train{1}, 1);
% ↓ two hidden layers
numHiddenUnits1 = numFeatures*6; % try this as OptParameter
numHiddenUnits2 = numFeatures*3; % try this as OptParameter
% ↓ feature dimension of each frame of output data
numResponses = size(y_train{1}, 1);
% ↓ para initialization
% dropoutLayer1=0.2; % try this as OptParameter
dropout_rate = 0.15; % try this as OptParameter
NEpoch=80; % no need to optimise, simply allow for a generous max
LearningRate=0.002; % try this as OptParameter (useful but secondary. Worth tuning within a narrow log range (1e-5 … 1e-3)).
% shuffle = 'once';
shuffle = 'every-epoch'; % sequence-to-sequence training and estimation
% layer1='last';
layer1='sequence';
mini_batch_size=128; % try this as OptParameter (lower priority, try options {64, 128, 256})
L2Regularization = 1e-04; % try this as OptParameter
GradientThreshold = 1;  % Rarely improves results when tuned, unless exploding gradients are a clear issue. Leave fixed (1 or 2).
LearnRateDropPeriod = 4; % % These fine-tune convergence schedule, but Bayesian optimisation on them is inefficient because their effect is weaker than LR/units/dropout.
LearnRateDropFactor = 0.6; 
ValidationPatience = 12;

% -------------------------
% Bayesian optimisation setup
% -------------------------
numHiddenUnits1_array = [round(0.25*numFeatures), round(0.5*numFeatures), numFeatures, 2*numFeatures, 4*numFeatures, 8*numFeatures];
numHiddenUnits2_array = [round(0.25*numFeatures), round(0.5*numFeatures), numFeatures, 2*numFeatures, 4*numFeatures, 8*numFeatures];
dropout_rate_array = [0.10, 0.15, 0.2, 0.25, 0.30]*100;
LearningRate_array = [1e-05, 5e-05, 1e-4, 5e-04, 1e-03, 5e-03, 1e-2]*1e05;
L2Regularization_array =[1e-5, 1e-04, 1e-3]*1e5;
mini_batch_size_array = [64, 128, 256];
penaliseValue = (10-allMean_tc_l)/allStd_tc_l;

optimArrays= struct('numHiddenUnits1',numHiddenUnits1_array,'numHiddenUnits2',numHiddenUnits2_array,'dropout_rate',dropout_rate_array,...
    'LearningRate',LearningRate_array,'L2Regularization',L2Regularization_array,'mini_batch_size',mini_batch_size_array,'penaliseValue',penaliseValue);

optimVars = [
    optimizableVariable('numHiddenUnits1_idx', [1, numel(numHiddenUnits1_array)], 'Type', 'integer')
    optimizableVariable('numHiddenUnits2_idx', [1, numel(numHiddenUnits2_array)], 'Type', 'integer')
    optimizableVariable('dropout_rate_idx',    [1, numel(dropout_rate_array)], 'Type', 'integer')
    optimizableVariable('LearningRate_idx',    [1, numel(LearningRate_array)], 'Type', 'integer')
    optimizableVariable('L2Regularization_idx',[1, numel(L2Regularization_array)], 'Type', 'integer')
    optimizableVariable('mini_batch_size_idx', [1, numel(mini_batch_size_array)], 'Type', 'integer')
];

% Objective: minimise validation RMSE (computed from predictions on x_val,y_val)
ObjFcn = @(T) objectiveFcn_seq2seq(T, x_train, y_train, x_val, y_val, ...
                                   numFeatures, numResponses, NEpoch, ValidationPatience,optimArrays);
                               
MaxEvals  = 50;   % adjust up if you have more compute

results = bayesopt(ObjFcn, optimVars, ...
    'MaxObjectiveEvaluations', MaxEvals, ...
    'IsObjectiveDeterministic', false, ...
    'UseParallel', false, ...      % set true if you have a pool/GPU capacity
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, ...
    'AcquisitionFunctionName', 'expected-improvement-plus');

best = results.XAtMinObjective;
disp('Best hyperparameters found by BO:');
disp(best);

%% Train the model on the best parameters obtained from Bayesian
% optimisation
[layersBest, optionsBest] = buildModelAndOptions(best, numFeatures, numResponses, ...
    NEpoch, ValidationPatience, x_test, y_test,optimArrays);

[net, info] = trainNetwork(x_train, y_train, layersBest, optionsBest);

% Evaluate on validation & test
valRMSE  = rmse_seqcell(predict(net, x_val), y_val);
testRMSE = rmse_seqcell(predict(net, x_test), y_test);

fprintf('Validation RMSE (best): %.4f\n', valRMSE);
fprintf('Test RMSE (best):       %.4f\n', testRMSE);

%% calculate RMSE, R2, Bias, 95%CI

% scale output back
% y_true_train = cellfun(@(v) v ./ scale_y, y_train, 'UniformOutput', false);
% y_pred_train = cellfun(@(v) v ./ scale_y, predict(net, x_train), 'UniformOutput', false);
% y_true_val = cellfun(@(v) v ./ scale_y, y_val, 'UniformOutput', false);
% y_pred_val = cellfun(@(v) v ./ scale_y, predict(net, x_val), 'UniformOutput', false);

y_true_train = cellfun(@(v) v .* allStd_tc_l + allMean_tc_l, y_train, 'UniformOutput', false);
y_pred_train = cellfun(@(v) v .* allStd_tc_l + allMean_tc_l, predict(net, x_train), 'UniformOutput', false);
y_true_val   = cellfun(@(v) v .* allStd_tc_l + allMean_tc_l, y_val, 'UniformOutput', false);
y_pred_val   = cellfun(@(v) v .* allStd_tc_l + allMean_tc_l, predict(net, x_val), 'UniformOutput', false);

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
fprintf('Train R2 = %.4f \n', R2_train);
fprintf('Train Bias = %.2f mm\n', bias_train);
fprintf('Train 95CI = [%.2f mm, %.2f mm]\n', LoA_lower_train, LoA_upper_train);

% [VALIDATION DATASET] calculate RMSE
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
fprintf('Validation R2 = %.4f \n', R2_val);
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

%% Functions

function [xOut, yOut] = dropNaNCells(xIn, yIn)
    emptyIdx = false(numel(yIn), 1);
    for k = 1:numel(yIn)
        yi = yIn{k};
        xi = xIn{k};
        if isempty(yi) || isempty(xi)
            emptyIdx(k) = true; continue
        end
        validIdx = ~any(isnan(yi), 1);
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

function [layers, options] = buildModelAndOptions(T, numFeatures, numResponses, NEpoch, ValidationPatience, x_test, y_test,optimArrays)

% Static Parameters
GradientThreshold = 1;  % Rarely improves results when tuned, unless exploding gradients are a clear issue. Leave fixed (1 or 2).
LearnRateDropPeriod = 4; % % These fine-tune convergence schedule, but Bayesian optimisation on them is inefficient because their effect is weaker than LR/units/dropout.
LearnRateDropFactor = 0.6; 
shuffle = 'every-epoch'; % sequence-to-sequence training and estimation

% Opt parameters
numHiddenUnits1 = optimArrays.numHiddenUnits1(T.numHiddenUnits1_idx);
numHiddenUnits2 = optimArrays.numHiddenUnits2(T.numHiddenUnits2_idx);
dropout_rate = optimArrays.dropout_rate(T.dropout_rate_idx)/100;
LearningRate = optimArrays.LearningRate(T.LearningRate_idx)*1e-05;
L2Regularization = optimArrays.L2Regularization(T.L2Regularization_idx)*1e-05;
mini_batch_size = optimArrays.mini_batch_size(T.mini_batch_size_idx);

layers = [
    sequenceInputLayer(numFeatures) % input layer
    bilstmLayer(numHiddenUnits1, 'OutputMode', 'sequence')  % 1st layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(dropout_rate)

    lstmLayer(numHiddenUnits2, 'OutputMode', 'sequence')  % 2nd layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(dropout_rate)

    fullyConnectedLayer(numResponses) % output layer
    regressionLayer
    ];

options = trainingOptions('adam',...
    'ValidationData',{x_test, y_test},...
    'ValidationPatience', ValidationPatience,... % changed
    'MiniBatchSize', mini_batch_size,...
    'MaxEpochs', NEpoch, ...
    'InitialLearnRate', LearningRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', LearnRateDropPeriod, ... % changed
    'LearnRateDropFactor', LearnRateDropFactor, ... % changed
    'GradientThreshold', GradientThreshold, ... % added
    'L2Regularization', L2Regularization, ...
    'Shuffle', shuffle, ...
    'Verbose', false, ...
    'Plots', 'training-progress',...
    'ExecutionEnvironment', 'gpu');
end


function results = objectiveFcn_seq2seq(T, x_train, y_train, x_test, y_test, ...
                                        numFeatures, numResponses, NEpoch, ValidationPatience,optimArrays)
    % Build model for this hyperparameter sample
    [layers, options] = buildModelAndOptions(T, numFeatures, numResponses, NEpoch, ValidationPatience, x_test, y_test,optimArrays);

    % Train
    try
        net = trainNetwork(x_train, y_train, layers, options);
    catch ME
        warning(ME.identifier, '%s', ME.message);
        results = optimArrays.penaliseValue; % penalize failure with '10' which is not too large like inf and therefor will not mess up gradients too much
%         results = table(optimArrays.penaliseValue, 'VariableNames', {'Objective'}); % penalize failure with '10' which is not too large like inf and therefor will not mess up gradients too much
        return
    end

    % Validate: compute RMSE on x_val,y_val
    y_pred_val = predict(net, x_test, 'MiniBatchSize', options.MiniBatchSize);
    results    = rmse_seqcell(y_pred_val, y_test);

    % Return objective (minimise RMSE)
%     results = table(rmseVal, 'VariableNames', {'Objective'});
end


function e = rmse_seqcell(YhatCell, YtrueCell)
    % Concatenate over sequences & timesteps; supports sequence-to-sequence (1 x T) or (D x T)
    assert(numel(YhatCell)==numel(YtrueCell), 'Prediction/target cell sizes must match.');
    seSum = 0; nTot = 0;
    for i = 1:numel(YhatCell)
        Yh = YhatCell{i};
        Yt = YtrueCell{i};
        % Ensure shapes match
        if ~isequal(size(Yh), size(Yt))
            % try to align (e.g., column vectors vs row vectors issues)
            if isequal(size(Yh,1), size(Yt,1)) && isequal(size(Yh,2), size(Yt,2))
                % ok
            else
                error('Size mismatch in sequence %d: pred %s vs true %s', ...
                    i, mat2str(size(Yh)), mat2str(size(Yt)));
            end
        end
        seSum = seSum + sum((Yh - Yt).^2, 'all');
        nTot  = nTot  + numel(Yh);
    end
    e = sqrt(seSum / max(1,nTot));
end


%% Save figs and results
% Define results path
local_data_dir=fullfile(pwd, 'TrainingModelsResults');
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

%% Save figures and useful script information
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

ModelTuningParams = struct('NEpoch',NEpoch,'numFeatures',numFeatures,'numResponses',numResponses,...
    'window',window,'step_between_points',step_between_points,'x_set',x_set);
ResultsParams = struct('model_net',net,'model_info',info,'bias_train',bias_train,'bias_val',bias_val,'LoA_lower_train',LoA_lower_train,'LoA_upper_train',LoA_upper_train,'LoA_lower_val',LoA_lower_val,'LoA_upper_val',LoA_upper_val,'R2_train',R2_train,'R2_val',R2_val,'train_rmse',train_rmse,'val_rmse',val_rmse);

save([local_data_dir filesep 'Trial_details.mat'],'ModelTuningParams','ResultsParams')
