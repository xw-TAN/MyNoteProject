window = 60; % approximately one cycle (data window length)
step_between_points = round(window/10);
averageMaxWindowDuration = 1.18; % [s]

% to keep consistency [comment these since involving multi-stage training]
rng(123, 'threefry'); % for randperm
gpurng(123, 'threefry'); % for model training

load("postZscore_withoutTimeNormalization_03_Sep_onlyEdgeCase.mat");

% proportion of edge data
edgeProportion = 0.20; % note: size ratio of all edge to all normal is 1:100

% percentage of datasets for different purposes
percentTra = 0.7;
percentVal = 0.2;
percentTes = 0.1;

%% Edge Case Data
x_total = [];
y_total = [];
totalLen= length(nAll_OutputE.tc);
for i = 1:step_between_points:totalLen-window
    % show the percent progress
    index = i:i+window-1;
    disp(['edge progress: ' num2str(i/(totalLen-window)*100, '%.3f%%')]);

    % output
    y_set = {nAll_OutputE.tc(index)'};

    % input
    temp = nAll_InputE.shoeTime(index, :);
    if max(temp) - min(temp) > 5 % in case of the window crossing two trials
        continue;
    end
    temp = (temp - temp(1)) / averageMaxWindowDuration;

    x_set = {[
        temp'; ...
        nAll_InputE.accel(index, :)'; ...
        nAll_InputE.EulerAngles_footIMU(index, [1,3])';...
        nAll_InputE.gyro(index, :)';...
        nAll_InputE.inc(index, :)';...
        nAll_InputE.shoeSize(index, :)';...
        nAll_InputE.thumb(index, :)';...
        nAll_InputE.legSide(index, :)']};

    x_total = [x_total; x_set];
    y_total = [y_total; y_set];
end

[x_total, y_total] = dropNaNCells(x_total, y_total);

lengthCell = length(y_total);
length_tra = floor(lengthCell * percentTra);
length_val = floor(lengthCell * percentVal);
length_tes = floor(lengthCell * percentTes);

% shuffle
idx = randperm(lengthCell);
x_total = x_total(idx);
y_total = y_total(idx);

% separate datasets
x_traE = x_total(1:length_tra);
y_traE = y_total(1:length_tra);
x_valE = x_total(length_tra+1:length_tra+length_val);
y_valE = y_total(length_tra+1:length_tra+length_val);
x_tesE = x_total(end-length_tes+1:end);
y_tesE = y_total(end-length_tes+1:end);


%% Normal Data
% x_total = [];
% y_total = [];
% totalLen= length(nAll_Output.tc);
% for i = 1:step_between_points:totalLen-window
%     % show the percent progress
%     index = i:i+window-1;
%     disp(['normal progress: ' num2str(i/(totalLen-window)*100, '%.3f%%')]);
% 
%     % output
%     y_set = {nAll_Output.tc(index)'};
% 
%     % input
%     temp = nAll_Input.shoeTime(index, :);
%     if max(temp) - min(temp) > 5 % in case of the window crossing two trials
%         continue;
%     end
%     temp = (temp - temp(1)) / averageMaxWindowDuration;
% 
%     x_set = {[
%         temp'; ...
%         nAll_Input.accel(index, :)'; ...
%         nAll_Input.EulerAngles_footIMU(index, [1,3])';...
%         nAll_Input.gyro(index, :)';...
%         nAll_Input.inc(index, :)';...
%         nAll_Input.shoeSize(index, :)';...
%         nAll_Input.thumb(index, :)';...
%         nAll_Input.legSide(index, :)']};
% 
%     x_total = [x_total; x_set];
%     y_total = [y_total; y_set];
% end
% 
% [x_total, y_total] = dropNaNCells(x_total, y_total);

load("Normal_xTotal_yTotal_shuffled_removedTestData.mat");

lengthCell = length(y_total);
length_tra = floor(lengthCell * percentTra);
length_val = floor(lengthCell * percentVal);
length_tes = floor(lengthCell * percentTes);

% shuffle
idx = randperm(lengthCell);
x_total = x_total(idx);
y_total = y_total(idx);

% separate datasets
x_traN = x_total(1:length_tra);
y_traN = y_total(1:length_tra);
x_valN = x_total(length_tra+1:length_tra+length_val);
y_valN = y_total(length_tra+1:length_tra+length_val);
x_tesN = x_total(end-length_tes+1:end);
y_tesN = y_total(end-length_tes+1:end);


%% Generate edge-normal mixed dataset
lengthN = length([x_traN; x_valN; x_tesN]);
lengthE = length([x_traE; x_valE; x_tesE]);

intendLengthN = lengthE / edgeProportion * (1-edgeProportion);
intendLengthTra = floor(length(x_traN) * intendLengthN/lengthN);
intendLengthVal = floor(length(x_valN) * intendLengthN/lengthN);
intendLengthTes = floor(length(x_tesN) * intendLengthN/lengthN);

x_traN = x_traN(1:intendLengthTra);
y_traN = y_traN(1:intendLengthTra);
x_valN = x_valN(1:intendLengthVal);
y_valN = y_valN(1:intendLengthVal);
x_tesN = x_tesN(1:intendLengthTes);
y_tesN = y_tesN(1:intendLengthTes);

% merge
x_tra =[x_traN; x_traE];
y_tra =[y_traN; y_traE];
x_val =[x_valN; x_valE];
y_val =[y_valN; y_valE];
x_tes =[x_tesN; x_tesE];
y_tes =[y_tesN; y_tesE];

% shuffle
idx_tra = randperm(length(x_tra));
x_tra = x_tra(idx_tra);
y_tra = y_tra(idx_tra);

idx_val = randperm(length(x_val));
x_val = x_val(idx_val);
y_val = y_val(idx_val);

idx_tes = randperm(length(x_tes));
x_tes = x_tes(idx_tes);
y_tes = y_tes(idx_tes);

%% LSTMs
NEpoch_ft = 20;             % no need to have too many epoch
LearningRate_ft = 1e-3;     % ten times smaller than that in normal training
mini_batch_size_ft = 256;

options = trainingOptions('adam', ...
    'ValidationData',{x_val, y_val}, ...
    'ValidationPatience', 8, ...
    'MiniBatchSize', mini_batch_size_ft, ...
    'MaxEpochs', NEpoch_ft, ...
    'InitialLearnRate', LearningRate_ft, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ...
    'LearnRateDropFactor', 0.6, ...
    'GradientThreshold', 1, ...
    'L2Regularization', 1e-3, ...
    'Shuffle','every-epoch', ...
    'Verbose',false, ...
    'Plots','training-progress', ...
    'ExecutionEnvironment','auto', ...
    'OutputNetwork','best-validation-loss');  % retain the best one

load("net_base.mat"); % trained using normal data
if ~exist('net', 'var')
    layers_ft = net_base.Layers;
else
    net_last = net; % trained from last stage
    layers_ft = net_last.Layers;
end
[net, info] = trainNetwork(x_tra, y_tra, layers_ft, options);


%% calculate RMSE, R2, Bias, 95%CI
% restore the sequence
[~, idx_tra_sort] = sort(idx_tra);
x_tra = x_tra(idx_tra_sort);
y_tra = y_tra(idx_tra_sort);

[~, idx_tes_sort] = sort(idx_tes);
x_tes = x_tes(idx_tes_sort);
y_tes = y_tes(idx_tes_sort);

% fprintf('\nnet_base normal train dataset:\n');
% [y_true_netBase_traN, y_pred_netBase_traN] = ...
%     printIndicator(x_traN, y_traN, allStd_tc, allMean_tc, net_base);

% fprintf('\nnet_base normal test dataset:\n');
% [y_true_netBase_tesN, y_pred_netBase_tesN] = ...
%     printIndicator(x_tesN, y_tesN, allStd_tc, allMean_tc, net_base);

% fprintf('\nnet_base edge train dataset:\n');
% [y_true_netBase_traE, y_pred_netBase_traE] = ...
%     printIndicator(x_traE, y_traE, allStd_tc, allMean_tc, net_base);

% fprintf('\nnet_base edge test dataset:\n');
% [y_true_netBase_tesE, y_pred_netBase_tesE] = ...
%     printIndicator(x_tesE, y_tesE, allStd_tc, allMean_tc, net_base);

% fprintf('\nnet normal train dataset:\n');
% [y_true_net_traN, y_pred_net_traN] = ...
%     printIndicator(x_traN, y_traN, allStd_tc, allMean_tc, net);

fprintf('\nnet normal test dataset:\n');
[y_true_net_tesN, y_pred_net_tesN] = ...
    printIndicator(x_tesN, y_tesN, allStd_tc, allMean_tc, net);

% fprintf('\nnet edge train dataset:\n');
% [y_true_net_traE, y_pred_net_traE] = ...
%     printIndicator(x_traE, y_traE, allStd_tc, allMean_tc, net);

fprintf('\nnet edge test dataset:\n');
[y_true_net_tesE, y_pred_net_tesE] = ...
    printIndicator(x_tesE, y_tesE, allStd_tc, allMean_tc, net);


%% plot figures
% plotFigures(y_true_netBase_tesN, y_pred_netBase_tesN, window, step_between_points, 'netBase testN dataset');
% plotFigures(y_true_netBase_tesE, y_pred_netBase_tesE, window, step_between_points, 'netBase testE dataset');
% plotFigures(y_true_net_tesN, y_pred_net_tesN, window, step_between_points, 'net testN dataset');
% plotFigures(y_true_net_tesE, y_pred_net_tesE, window, step_between_points, 'net testE dataset');


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

function [y_true, y_pred] = printIndicator(x_set, y_set, allStd_tc, allMean_tc, net)
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