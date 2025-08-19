% Model Training

percentVal = [0.0, 0.1]; % 0-10% for validation
percentTra = [0.1, 0.8]; % 10-80% for training
percentTes = [0.8, 1.0]; % 80-100% for testing

window=60; % approximately one cycle (data window length)

step_between_points=window/10;   % sliding window
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



% codes below are unchanged except removing the y_scale factor
% (maybe codes below can be simplified using the above numerical vectors 
% instead of cell vectors)

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
% scale_y = 0.1;
% scale_y = 0.02;
% y_train = cellfun(@(v) v * scale_y, y_train, 'UniformOutput', false);
% y_test  = cellfun(@(v) v * scale_y, y_test,  'UniformOutput', false);
% y_val   = cellfun(@(v) v * scale_y, y_val,   'UniformOutput', false);

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

layers = [
    sequenceInputLayer(numFeatures) % input layer
    bilstmLayer(numFeatures*6, 'OutputMode', 'sequence')  % 1st layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(0.15)

    lstmLayer(numFeatures*3, 'OutputMode', 'sequence')  % 2nd layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(0.15)

    fullyConnectedLayer(numResponses) % output layer
    regressionLayer
    ];

options = trainingOptions('adam',...
    'ValidationData',{x_test, y_test},...
    'ValidationPatience', 12,... % changed
    'MiniBatchSize', mini_batch_size,...
    'MaxEpochs', NEpoch, ...
    'InitialLearnRate', LearningRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ... % changed
    'LearnRateDropFactor', 0.6, ... % changed
    'GradientThreshold', 1, ... % added
    'L2Regularization', 1e-4, ...
    'Shuffle', shuffle, ...
    'Verbose', false, ...
    'Plots', 'training-progress',...
    'ExecutionEnvironment', 'auto');


[net, info] = trainNetwork(x_train, y_train, layers, options);




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

% %% Save figs and results
% % Define results path
% local_data_dir=fullfile(pwd, 'TrainingModelsResults');
% % Observe present folders (named in arithmetic order)
% files = dir(local_data_dir);
% dirFlags= [files.isdir];
% subFolder = files(dirFlags);
% subFolderNames = [str2double({subFolder(3:end).name})];
% % Observe latest folder of results (largest number)
% maxFolderNum=max(subFolderNames);
% % Create the next folder
% nextFolderName = maxFolderNum+1;
% if nextFolderName<10
%     nextFolderName=['00' num2str(nextFolderName)];
% elseif nextFolderName<100
%     nextFolderName=['0' num2str(nextFolderName)];
% else
%     nextFolderName=['0' num2str(nextFolderName)];
% end
% mkdir(fullfile(local_data_dir,nextFolderName))
% local_data_dir=fullfile(local_data_dir,nextFolderName);
% 
% % Save figures and useful script information
% filename='TrainCorrFig.png';
% saveas(train_corr_fig,[local_data_dir filesep filename])
% filename='TrainCorrFig.fig';
% saveas(train_corr_fig,[local_data_dir filesep filename])
% 
% filename='TrainCompFig.png';
% saveas(train_comp_fig,[local_data_dir filesep filename])
% filename='TrainCompFig.fig';
% saveas(train_comp_fig,[local_data_dir filesep filename])
% 
% filename='ValCorrFig.png';
% saveas(val_corr_fig,[local_data_dir filesep filename])
% filename='ValCorrFig.fig';
% saveas(val_corr_fig,[local_data_dir filesep filename])
% 
% filename='ValCompFig.png';
% saveas(val_comp_fig,[local_data_dir filesep filename])
% filename='ValCompFig.fig';
% saveas(val_comp_fig,[local_data_dir filesep filename])
% 
% ModelTuningParams = struct('dropoutLayer1',dropoutLayer1,'LearningRate',LearningRate,'mini_batch_size',mini_batch_size,'NEpoch',NEpoch,'numFeatures',numFeatures,'numHiddenUnits1',numHiddenUnits1,'numHiddenUnits2',numHiddenUnits2,'numResponses',numResponses,...
%     'layers',layers,'options',options,'scale_y',scale_y,'window',window,'step_between_points',step_between_points,'x_set',x_set);
% ResultsParams = struct('model_net',net,'model_info',info,'bias_train',bias_train,'bias_val',bias_val,'LoA_lower_train',LoA_lower_train,'LoA_upper_train',LoA_upper_train,'LoA_lower_val',LoA_lower_val,'LoA_upper_val',LoA_upper_val,'R2_train',R2_train,'R2_val',R2_val,'train_rmse',train_rmse,'train_test_subjects',train_test_subjects,'val_rmse',val_rmse,'val_subjects',val_subjects);
% 
% save([local_data_dir filesep 'Trial_details.mat'],'ModelTuningParams','ResultsParams')
