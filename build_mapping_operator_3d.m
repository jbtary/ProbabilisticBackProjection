% Code to build the mapping operator to locate events: 3D + constant velocity
clear
addpath('prob_back_loc')

% Parameters needed
allsta = 'OBS_positions.txt'; % File with all station names and locations
vel = 1.5; % Velocity in km/s
rang = [0 50 0.2;0 50 0.2;0 3.2 0.2]; % Search range -- min, max, increments in km
grid0 = [3.0436111 -84]; % Define the 0,0 of the grid (lower-left corner)
[minx,miny] = deg2utm(grid0(1),grid0(2));
minx = minx/1000; miny = miny/1000;
Fs = 100; % Sampling frequency in Hz
lag = 1900; % Maximum cross-corr lag to consider in samples
tlag1 = -lag:lag;

% Name mapping operator to save
name_mapping = 'mapping_operator_whales_3d';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating maximum differential traveltimes for each pair')
% Get info for all stations
stations = importdata(allsta,' ');
sta_name = stations.textdata;
[sta(:,1),sta(:,2)] = deg2utm(stations.data(:,1),stations.data(:,2)); % lat lon
sta = sta/1000; % in km
sta(:,1) = sta(:,1) - minx; % Set at same 0,0 (not mandatory)
sta(:,2) = sta(:,2) - miny;
sta(:,3) = -stations.data(:,3); % Depth positive
Nsta = length(sta_name);

[ pair1, Npair1 ] = cons_pair1( Nsta ); % Get pairs of stations

for ii = 1:length(pair1(:,1))
    % Use constant vel
    dist_diff = sqrt( ( sta(pair1(ii,1),1) - sta(pair1(ii,2),1) )^2 + ...
        ( sta(pair1(ii,1),2) - sta(pair1(ii,2),2) )^2 + ...
        ( sta(pair1(ii,1),3) - sta(pair1(ii,2),3) )^2);
    max_tlag(ii) = dist_diff/vel;
    clear dist_diff
end
max_samlag = round(max_tlag*Fs);

% -------------------------------------------------------------------------
% Compute the no. of grid points in each direction
% -------------------------------------------------------------------------
x_no = (rang(1,2) - rang(1,1))/rang(1,3) + 1;
y_no = (rang(2,2) - rang(2,1))/rang(2,3) + 1;
z_no = (rang(3,2) - rang(3,1))/rang(3,3) + 1;

% -------------------------------------------------------------------------
% Compute the mapping operator for each pair
% -------------------------------------------------------------------------
disp('Computing backpropagation maps, slow...')
[CC1_ind,alpha] = mapping_var_3d(x_no,y_no,z_no,pair1,Npair1,...
    rang(1,1):rang(1,3):rang(1,2),rang(2,1):rang(2,3):rang(2,2),rang(3,1):rang(3,3):rang(3,2),...
    sta,vel,Fs,tlag1);

% Save the mapping operator for future uses ('-v7.3' because files get LARGE)
%alpha = single(alpha); CC1_ind = single(CC1_ind); % To save space 
save(name_mapping,'max_samlag','x_no','y_no','z_no','CC1_ind','alpha',...
        'rang','pair1','Npair1','sta','sta_name','tlag1','Fs','vel','-v7.3')
