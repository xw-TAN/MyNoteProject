function [euler_angles] = calc3DShankOrientation(s_upp, s_med, s_lat, side)
% input:
%   s_upp: shank_upper, Nx3，N frames and 3D coodinates
%   s_med: shank_med, Nx3
%   s_lat: shank_lat, Nx3
%   side: 1 char, 'L' or 'R', naming the leg side the markers placed on
%
% output：
%   euler_angles: Nx3，N frames (deg)
%

N = size(s_upp, 1);
R_all = zeros(3,3,N);

%% create marker-based frame [build a frame similar to Vicon frame]
for i = 1:N
    if side == 'L'
        % z-axis, right direction whether the left or right leg
        Z = s_med(i,:) - s_lat(i,:); 
        Z = Z / norm(Z);
    elseif side == 'R'
        Z = s_lat(i,:) - s_med(i,:);
        Z = Z / norm(Z);
    else
        disp('side value input error');
    end
   
    % midpoint
    mid_front = (s_med(i,:) + s_lat(i,:)) / 2;

    temp = cross(s_upp(i,:) - mid_front, Z); % x-axis, forward direction
    X = temp / norm(temp);
    
    % y-axis, upward direction
    Y = cross(Z, X);    
    Y = Y / norm(Y);

    % rotationay matrix of marker-based frame relative to Vicon world frame
    % the code below ensures it represents the relation of Marker to Vicon
    R_all(:,:,i) = [X' Y' Z'];
end


%% calculate angles of marker-based frame to the Vicon frame

euler_angles = zeros(N,3);
for i = 1:N
    % representation of Marker frame relative to Vicon frame means to 
    % rotate Vicon frame to get Marker frame (From_Vicon_to_Marker)
    R_From_Vicon_to_Marker = R_all(:,:,i);
    
    % convert to euler angles (deg). X is least likely to cause dead-lock.
    % Z: flex/ext, X: inversion/eversion, Y: internal/external rotation
    euler_angles(i,:) = rad2deg(rotm2eul(R_From_Vicon_to_Marker, 'ZXY'));
end

% avoid angle jumps at +-180deg
euler_angles = unwrap(deg2rad(euler_angles), [], 1) * 180/pi;

end

