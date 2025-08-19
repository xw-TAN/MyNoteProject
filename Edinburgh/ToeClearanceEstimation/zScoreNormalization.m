total_subject = 1:10;

All_Output.tc_l = [];
All_Input.EulerAngles_footIMU_l = [];
All_Input.shoeTime = [];
All_Input.thumb_l  = [];
All_Input.little_l = [];
All_Input.heel_l   = [];
All_Input.accel_l  = [];
All_Input.gyro_l   = [];
All_Input.inc      = [];
All_Input.shoeSize = [];

for i = total_subject
    subjects = {'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'S10'};
    subject_shoe_size   = [42, 40, 42, 42, 42, 40, 40, 42, 40, 42]; % [EU Size]
    load butterworth_3rd_FS100Hz_FC15Hz.mat; % butterworth filter, 3rd order, fs 100Hz, fc 15Hz

    steps_to_skip_start =300; % neglect the first 300 points after the first HS
    steps_to_skip_end   =5; % neglect the last 5 points
    
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

    experiments =fieldnames(results.(subjects{i}).ShokacData);

    % ↓ loop each condition
    for j=1:length(experiments)

        experiment_name=experiments{j};
        inc=sscanf(strrep(experiment_name(strfind(experiment_name,'n')+1:end),'_','-'), '%f');

        % skip edge case data
        if contains(experiment_name, 'edge')
            continue;
        end

        % ↓ get all data of this condition under 'ShokacData'
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

        % remove the row from index if there are wrong elements in this row
        sub_data = experiment_data.masked_data(index, 1:5);
        rows_to_keep = all(sub_data == 1, 2);
        index = index(rows_to_keep);

        % get the time using the latest index
        eq_vicon_time = experiment_data.time(index) - experiment_data.time(1) - experiment_data.time_diff;

        % output: toe clearance
        tc = interp1(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).Timesteps,...
                     results.(subjects{i}).MarkerData.([experiment_name '_Processed']).toe_clearance_l,...
            eq_vicon_time);
        All_Output.tc_l = [All_Output.tc_l; filtfilt(SOS, G, tc)];

        % input: EulerAngles_footIMU_l
        EulerAngles_footIMU_l = ...
            interp1(results.(subjects{i}).MarkerData.([experiment_name '_Processed']).Timesteps,...
                    results.(subjects{i}).MarkerData.([experiment_name '_Processed']).EulerAngles_footIMU_l',...
            eq_vicon_time);
        All_Input.EulerAngles_footIMU_l = [All_Input.EulerAngles_footIMU_l; EulerAngles_footIMU_l];

        % input: time
        time = experiment_data.time-experiment_data.time(1);
        All_Input.shoeTime = [All_Input.shoeTime; time(index)];

        % input: force
        All_Input.thumb_l  = [All_Input.thumb_l; experiment_data.thumb_l(index,  2:7)];
        All_Input.little_l = [All_Input.little_l; experiment_data.little_l(index,  2:7)];
        All_Input.heel_l   = [All_Input.heel_l; experiment_data.heel_l(index,  2:7)];

        % input: acc
        All_Input.accel_l  = [All_Input.accel_l; experiment_data.accel_l(index, 1:3)];

        % input: gyro
        All_Input.gyro_l   = [All_Input.gyro_l; experiment_data.gyro_l(index, 1)];

        % input: inc
        All_Input.inc      = [All_Input.inc; inc*ones( length(experiment_data.gyro_l(index, 1)), 1)];

        % input: shoe size
        All_Input.shoeSize = [All_Input.shoeSize; subject_shoe_size(i)*ones( length(experiment_data.gyro_l(index, 1)), 1)];

        % display the current progress
        disp([char(subjects{i}) ': ' char(experiments{j})]);
    end
end

% calculate mean of all types of input/output data
allMean_tc_l = mean(All_Output.tc_l);
allMean_EulerAngles_footIMU_l = mean(All_Input.EulerAngles_footIMU_l);
allMean_shoeTime = mean(All_Input.shoeTime);
allMean_thumb_l  = mean(All_Input.thumb_l);
allMean_little_l = mean(All_Input.little_l);
allMean_heel_l   = mean(All_Input.heel_l);
allMean_accel_l  = mean(All_Input.accel_l);
allMean_gyro_l   = mean(All_Input.gyro_l);
allMean_inc      = mean(All_Input.inc);
allMean_shoeSize = mean(All_Input.shoeSize);

% calculate std of all types of input/output data
allStd_tc_l = std(All_Output.tc_l);
allStd_EulerAngles_footIMU_l = std(All_Input.EulerAngles_footIMU_l);
allStd_shoeTime = std(All_Input.shoeTime);
allStd_thumb_l  = std(All_Input.thumb_l);
allStd_little_l = std(All_Input.little_l);
allStd_heel_l   = std(All_Input.heel_l);
allStd_accel_l  = std(All_Input.accel_l);
allStd_gyro_l   = std(All_Input.gyro_l);
allStd_inc      = std(All_Input.inc);
allStd_shoeSize = std(All_Input.shoeSize);

% z-score normalization
nAll_Output.tc_l    = (All_Output.tc_l - allMean_tc_l)./allStd_tc_l;
nAll_Input.EulerAngles_footIMU_l = (All_Input.EulerAngles_footIMU_l - allMean_EulerAngles_footIMU_l)./allStd_EulerAngles_footIMU_l;
nAll_Input.shoeTime = (All_Input.shoeTime - allMean_shoeTime)./allStd_shoeTime;
nAll_Input.thumb_l  = (All_Input.thumb_l - allMean_thumb_l)./allStd_thumb_l;
nAll_Input.little_l = (All_Input.little_l - allMean_little_l)./allStd_little_l;
nAll_Input.heel_l   = (All_Input.heel_l - allMean_heel_l)./allStd_heel_l;
nAll_Input.accel_l  = (All_Input.accel_l - allMean_accel_l)./allStd_accel_l;
nAll_Input.gyro_l   = (All_Input.gyro_l - allMean_gyro_l)./allStd_gyro_l;
nAll_Input.inc      = (All_Input.inc - allMean_inc)./allStd_inc;
nAll_Input.shoeSize = (All_Input.shoeSize - allMean_shoeSize)./allStd_shoeSize;