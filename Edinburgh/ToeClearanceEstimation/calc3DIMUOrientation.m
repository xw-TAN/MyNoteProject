function [euler_angles] = calc3DIMUOrientation(heel, mtp1, mtp5, side)
% input:
%   heel: Nx3，N frames and 3D coodinates
%   mtp1: Nx3
%   mtp5: Nx3
%   side: 1 char, naming the leg side of markers
%
% output：
%   euler_angles: Nx3，N frames (deg)，X-Y-Z sequence
%

N = size(heel, 1);
R_all = zeros(3,3,N);

%% create marker-based frame
for i = 1:N
    if side == 'L'
        % x-axis from mtp1 to mtp5, pointing the left-side
        X = mtp5(i,:) - mtp1(i,:);
        X = X / norm(X);
    elseif side == 'R'
        X = mtp1(i,:) - mtp5(i,:); 
        X = X / norm(X);
    else
        disp('side value input error');
    end
   
    % midpoint between mtp1 and mtp5
    mid_front = (mtp1(i,:) + mtp5(i,:)) / 2;

    % cross vector between x-axis and the vector from heel to midpoint
    temp = cross(X, mid_front - heel(i,:));
    Z = temp / norm(temp); % downward

    % y-axis, pointing the forward direction
    Y = cross(Z, X);
    Y = Y / norm(Y);

    % rotationay matrix of marker-based frame relative to Vicon world frame
    R_all(:,:,i) = [X' Y' Z'];
end


%% calculate angles of IMU sensor frame in IMU world frame

% IMU world frame rotates -90deg around its x-axis to get the Vicon frame
R_From_IMUWorld_to_Vicon    = rotx(-90);
if side == 'L'
    R_From_Marker_to_IMUSensor  = roty(5) * rotx(5) * rotz(5); % right-multiplied (around self-frame)
elseif side == 'R'
    R_From_Marker_to_IMUSensor  = roty(-5) * rotx(5) * rotz(-5);
else
    disp('side value input error');
end

euler_angles = zeros(N,3);
for i = 1:N
    % representation of Marker frame relative to Vicon frame means to 
    % rotate Vicon frame to get Marker frame (From_Vicon_to_Marker)
    R_From_Vicon_to_Marker = R_all(:,:,i);
    
    % get the representation of IMU sensor frame relative to IMU world 
    % frame, so From_IMUWorld_to_IMUSensor
    R_From_IMUWorld_to_IMUSensor = R_From_IMUWorld_to_Vicon * R_From_Vicon_to_Marker * R_From_Marker_to_IMUSensor;

    % convert to euler angles in x-y-z sequency (deg). Y is least likely to
    % cause dead-lock
    euler_angles(i,:) = rad2deg(rotm2eul(R_From_IMUWorld_to_IMUSensor, 'XYZ'));
end

end

