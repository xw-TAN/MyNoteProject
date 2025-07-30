%
% To add new info into the old dataset (downloaded from Andreas oneDrive), 
% and form a new dataset for 'ReadingShoeDataTrainingModels_Xiaowei_v3.m'
%

clc;
clear;

% subject number
num = 10;

% loop each subject
for i = 10:num
    % ↓ load old dataset
    subNo = ['S' num2str(i)];
    load([subNo '.mat']);
    
    % ↓ wrap
    R.S = eval(subNo);
    
    fns = fieldnames( R.S.ShokacData );
    ConditionNum = numel(fns);
    
    % ↓ loop each condition
    for j = 1:ConditionNum

        % ↓ replace '_' with '-' to match the files' name
        fns_new = fns{j};
        underscores = find(fns_new == '_');
        if numel(underscores) == 2
            fns_new(underscores(2)) = '-';
        end



        % % ↓ get shoeData file path
        % pathName = ['./' subNo '/ShoeData' '/' fns_new '.txt'];
        % % ↓ get file data
        % shoeData = readmatrix(pathName);
        % % ↓ get time sequence
        % shoeTime = shoeData(:, 1);
        % % ↓ get left COM
        % shoeCOML = shoeData(:, 66:67);
        % % ↓ get right COM for HS identification
        % shoeCOMR = shoeData(:, 68:69);
        % % ↓ get the right HS index from shoeCOMR
        % idxHS_s = find( diff(sum(shoeCOMR, 2) ~= 0) == 1) + 1; % add 1 is for the latter one
        % ↓ add new info into the old dataset
        % R.S.ShokacData.(fns{j}).hs_r_time = shoeTime(idxHS_s);
        % R.S.ShokacData.(fns{j}).shokac_starting_time = shoeTime(idxHS_s(1)) - shoeTime(1);
        % R.S.ShokacData.(fns{j}).phase_l = sum(shoeCOML, 2) ~= 0; % 1 stance, 0 swing



        % ↓ get grfData file path
        % pathName = ['./' subNo '/GRFData' '/' fns_new '_GRFs_Processed.mot'];
        % % ↓ get file data
        % grfData = readmatrix(pathName, 'FileType', 'text');
        % % ↓ get the right vertical force for HS identification
        % grfData = grfData(:, [1, 3]);
        % % ↓ get the right HS time from Treadmill
        % idxHS_t = find( diff(grfData(:, 2) > 50) == 1) + 1; % add 1 is for the latter one
        % ↓ add new info into the old dataset
        % R.S.ShokacData.(fns{j}).time_diff = shoeTime(idxHS_s(1)) - shoeTime(1) - grfData(idxHS_t(1), 1);



        % ↓ get markerData file path
        pathName = ['./' subNo '/MarkerData' '/' fns_new '_Processed.trc'];
        % ↓ get file data
        markerData = readmatrix(pathName, 'FileType', 'text');
        % ↓ get the Timesteps
        markerTime = markerData(:, 2);
        % ↓ get the Left MTP1, MTP5, HEEL
        markerMTP1_L = markerData(:, 78:80);
        markerMTP5_L = markerData(:, 81:83);
        markerHEEL_L = markerData(:, 75:77);
        % ↓ get the Right MTP1, MTP5, HEEL
        markerMTP1_R = markerData(:, 66:68);
        markerMTP5_R = markerData(:, 69:71);
        markerHEEL_R = markerData(:, 63:65);
        % ↓ delete NaN elements
        idxNan = isnan(markerTime) | ...
                 any(isnan(markerMTP1_L), 2) | ...
                 any(isnan(markerMTP5_L), 2) | ...
                 any(isnan(markerHEEL_L), 2) | ...
                 any(isnan(markerMTP1_R), 2) | ...
                 any(isnan(markerMTP5_R), 2) | ...
                 any(isnan(markerHEEL_R), 2);
        markerTime(idxNan) = [];
        markerMTP1_L(idxNan, :) = [];
        markerMTP5_L(idxNan, :) = [];
        markerHEEL_L(idxNan, :) = [];
        markerMTP1_R(idxNan, :) = [];
        markerMTP5_R(idxNan, :) = [];
        markerHEEL_R(idxNan, :) = [];

        % ↓ add new info into the old dataset
        R.S.MarkerData.([fns{j} '_Processed']).data.Timesteps = markerTime';

        % ↓ calculate 3D IMU orientation angles
        euler_angles_l = calc3DIMUOrientation(markerHEEL_L, markerMTP1_L, markerMTP5_L, 'L');
        euler_angles_r = calc3DIMUOrientation(markerHEEL_R, markerMTP1_R, markerMTP5_R, 'R');
        % ↓ add new info into the old dataset
        R.S.MarkerData.([fns{j} '_Processed']).data.EulerAngles_l = euler_angles_l;
        R.S.MarkerData.([fns{j} '_Processed']).data.EulerAngles_r = euler_angles_r;


        % ↓ display the current progress
        disp([char(subNo) ': ' char(fns(j))]);
    end

    % override the old dataset using the new one
    eval([subNo ' = R.S']);
    vars = who('S*');
    clearvars('-except', vars{:});
end

