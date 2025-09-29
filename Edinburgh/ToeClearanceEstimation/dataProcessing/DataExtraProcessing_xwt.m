%
% To add new info into the old dataset and form a new dataset for use in
% 'ReadingShoeDataTrainingModels_Xiaowei_v3.m'
%

clc;
clear;

% subject number
num = 10;

% loop each subject
for i = 1:num
    local_data_dir=fullfile(pwd,'..','..','..','Moonshot','LocalData','Experiment','FESTrajectories_v02');
    
    % ↓ load old dataset
    subNo = ['S' num2str(i)];
    load([local_data_dir filesep subNo '.mat']);
    
    % ↓ wrap
    R.S = eval(subNo);

    % ↓ get condition numbers for this subject
    fns = fieldnames( R.S.ShokacData );
    ConditionNum = numel(fns);
    
    % ↓ loop each condition
    for j = 1:ConditionNum

        % ↓ replace '_' with '-' to match the files' name
        fns_new = fns{j}; % used to search the files
        underscores = find(fns_new == '_');
        if numel(underscores) == 2 && ~contains(fns_new,'edgecase')
            fns_new(underscores(2)) = '-';
        end


        % % ↓ get shoeData file path
        % pathName = [local_data_dir filesep subNo filesep 'ShoeData' filesep fns_new '.txt'];
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
        % pathName = [local_data_dir filesep subNo filesep 'GRFData' filesep fns_new '_GRFs_Processed.mot'];
        % % ↓ get file data
        % grfData = readmatrix(pathName, 'FileType', 'text');
        % % ↓ get the right vertical force for HS identification
        % grfData = grfData(:, [1, 3]);
        % % ↓ get the right HS time from Treadmill
        % idxHS_t = find( diff(grfData(:, 2) > 50) == 1) + 1; % add 1 is for the latter one
        % ↓ add new info into the old dataset
        % R.S.ShokacData.(fns{j}).time_diff = shoeTime(idxHS_s(1)) - shoeTime(1) - grfData(idxHS_t(1), 1);


        % ↓ get markerData file path
        if contains(fns_new,'edgecase')
            pathName = [local_data_dir filesep subNo filesep 'MarkerData' filesep 'EdgeCases' filesep fns_new '_Processed.trc'];
        else
            pathName = [local_data_dir filesep subNo filesep 'MarkerData' filesep fns_new '_Processed.trc'];
        end
        % ↓ get numerical data
        markerData = readmatrix(pathName, 'FileType', 'text');

        % delete rows with NaN
        rowsWithNaN = any(isnan(markerData), 2);
        markerData(rowsWithNaN, :) = [];
        
        % ↓ get title line
        lines = readlines(pathName);
        Title = split(lines(4));

        % ↓ check if there are duplicate marker names
        [uniqueTitle, ~, idx] = unique(Title);
        counts = histcounts(idx, 1:numel(uniqueTitle)+1);
        repeatedStrings = uniqueTitle(counts > 1);
        for n = 1:numel(uniqueTitle)
            if counts(n) > 1
                fprintf('%s: "%s" appears %d times\n', pathName, uniqueTitle{n}, counts(n));
            end
        end
        % if there are duplicate marker names, skip this file
        if length(uniqueTitle) ~= length(Title)
            continue;
        end
        
        % ↓ get the Timesteps
        markerTime = markerData(:, 2);

        % ↓ get the Left MTP1, MTP5, HEEL
        idx = find(contains(Title, 'L_Heel'));
        if isempty(idx)
            disp(['Error: [L_Heel] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerHEEL_L = markerData(:, rag);
        
        idx = find(contains(Title, 'L_MTP1'));
        if isempty(idx)
            disp(['Error: [L_MTP1] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerMTP1_L = markerData(:, rag);
        
        idx = find(contains(Title, 'L_MTP5'));
        if isempty(idx)
            disp(['Error: [L_MTP5] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerMTP5_L = markerData(:, rag);
        
        % ↓ get the Right MTP1, MTP5, HEEL
        idx = find(contains(Title, 'R_Heel'));
        if isempty(idx)
            disp(['Error: [R_Heel] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerHEEL_R = markerData(:, rag);
        
        idx = find(contains(Title, 'R_MTP1'));
        if isempty(idx)
            disp(['Error: [R_MTP1] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerMTP1_R = markerData(:, rag);
        
        idx = find(contains(Title, 'R_MTP5'));
        if isempty(idx)
            disp(['Error: [R_MTP5] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerMTP5_R = markerData(:, rag);

        % ↓ get the Left shank_upper, shank_med, shank_lat
        idx = find(contains(Title, 'L_Shank_Upper'));
        if isempty(idx)
            disp(['Error: [L_Shank_Upper] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerSKUP_L = markerData(:, rag);

        idx = find(contains(Title, 'L_Shank_Med'));
        if isempty(idx)
            disp(['Error: [L_Shank_Med] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerSKME_L = markerData(:, rag);

        idx = find(contains(Title, 'L_Shank_Lat'));
        if isempty(idx)
            disp(['Error: [L_Shank_Lat] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerSKLA_L = markerData(:, rag);

        % ↓ get the right shank_upper, shank_med, shank_lat
        idx = find(contains(Title, 'R_Shank_Upper'));
        if isempty(idx)
            disp(['Error: [R_Shank_Upper] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerSKUP_R = markerData(:, rag);

        idx = find(contains(Title, 'R_Shank_Med'));
        if isempty(idx)
            disp(['Error: [R_Shank_Med] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerSKME_R = markerData(:, rag);

        idx = find(contains(Title, 'R_Shank_Lat'));
        if isempty(idx)
            disp(['Error: [R_Shank_Lat] ' pathName]);
            continue;
        end
        rag = (idx-3)*3+3:(idx-3)*3+5;
        markerSKLA_R = markerData(:, rag);

        % ↓ add time steps
        R.S.MarkerData.([fns{j} '_Processed']).Timesteps = markerTime';

        % ↓ calculate 3D foot IMU orientation angles
        euler_angles_imu_l = calc3DIMUOrientation(markerHEEL_L, markerMTP1_L, markerMTP5_L, 'L');
        euler_angles_imu_r = calc3DIMUOrientation(markerHEEL_R, markerMTP1_R, markerMTP5_R, 'R');

        % ↓ calculate 3D shank orientation angles
        euler_angles_shank_l = calc3DShankOrientation(markerSKUP_L, markerSKME_L, markerSKLA_L, 'L');
        euler_angles_shank_r = calc3DShankOrientation(markerSKUP_R, markerSKME_R, markerSKLA_R, 'R');
        
        % ↓ add euler angles
        R.S.MarkerData.([fns{j} '_Processed']).EulerAngles_footIMU_l = euler_angles_imu_l';
        R.S.MarkerData.([fns{j} '_Processed']).EulerAngles_footIMU_r = euler_angles_imu_r';
        R.S.MarkerData.([fns{j} '_Processed']).EulerAngles_shank_l   = euler_angles_shank_l';
        R.S.MarkerData.([fns{j} '_Processed']).EulerAngles_shank_r   = euler_angles_shank_r';

        % ↓ display the current progress
        disp([char(subNo) ': ' char(fns(j))]);
    end

    % override the old dataset using the new one
    eval([subNo ' = R.S;']);
    vars = who('S*');
    clearvars('-except', vars{:});
end

