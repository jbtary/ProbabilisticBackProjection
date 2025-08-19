% Code to build the mapping operator to locate events using 3d travel-time tables
clear
addpath('prob_back_loc')

% Parameters needed
allsta = 'OBS_positions.txt'; % File with all station names and locations
path1 = 'TTimes/'; % Path to NonLinLoc travel-time (TT) tables
path2 = 'Velmod/'; % Path to corresponding NonLinLoc velocity model
wave = 'P'; % Wave type for TT tables and velocity model 
vpvs = 0; % >0 to use a Vp/Vs ratio instead of using S wave TT tables
rang = [0 50 0.2;0 50 0.2;0 3.2 0.2]; % Search range -- min, max, resolution in km
grid0 = [3.0436111 -84]; % Define the 0,0 of the grid (lower-left corner)
[minx,miny] = deg2utm(grid0(1),grid0(2));
minx = minx/1000; miny = miny/1000;
Fs = 100; % Sampling frequency in Hz
lag = 1900; % Maximum cross-corr lag to consider in samples
tlag1 = -lag:lag;

% Name mapping operator to save
name_mapping = ['mapping_operator_' wave '_whales'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating maximum differential traveltimes for each pair')
% Get info for all stations
stations = importdata(allsta,' ');
sta_name = stations.textdata;
[sta(:,1),sta(:,2)] = deg2utm(stations.data(:,1),stations.data(:,2)); % lat lon
sta = sta/1000; % in km
sta(:,1) = sta(:,1) - minx; % Set at same 0,0 (not mandatory)
sta(:,2) = sta(:,2) - miny;
sta(:,3) = stations.data(:,3);
if sum(stations.data(:,3)<0) > 0; sta(:,3) = -stations.data(:,3); end
Nsta = length(sta_name);

[ pair1, Npair1 ] = cons_pair1( Nsta ); % Get pairs of stations
[~,vel] = read_vel(path2,'velmod.P.mod');

for ii = 1:length(pair1(:,1))
    % Use travel-times from 3D Velocity model
    if ii > 1 && strcmp(sta_name{pair1(ii-1,1)},sta_name{pair1(ii,1)}) ~= 1
        clear ttsta1
        filename = ['ttimes.P.' sta_name{pair1(ii,1)} '.time'];
        [~,ttsta1,~,~] = read_tt(path1,filename);
    else
        filename = ['ttimes.P.' sta_name{pair1(ii,1)} '.time'];
        [~,ttsta1,~,~] = read_tt(path1,filename);
    end

    filename = ['ttimes.P.' sta_name{pair1(ii,2)} '.time'];
    [hdrgrid,ttsta2,~,~] = read_tt(path1,filename);
    if strcmp(wave,'S') == 1
        ttsta1 = ttsta1*vpvs; ttsta2 = ttsta2*vpvs;
    end
    max_tlag(ii) = diff_time_tt(ttsta1,ttsta2,hdrgrid,vel,rang(1,2),rang(2,2),rang(3,2));
    clear filename hdrgrid ttsta2
end
max_samlag = round(max_tlag*Fs); clear vel path2 ttsta1

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
[CC1_ind,alpha] = mapping_var_tt(x_no,y_no,z_no,pair1,Npair1,...
    rang(1,1):rang(1,3):rang(1,2),rang(2,1):rang(2,3):rang(2,2),rang(3,1):rang(3,3):rang(3,2),...
    sta,Fs,tlag1,path1,sta_name,wave,vpvs);

% Save the mapping operator for future uses ('-v7.3' because files get LARGE)
%alpha = single(alpha); CC1_ind = single(CC1_ind); % To save space 
save(name_mapping,'max_samlag','grid0','x_no','y_no','z_no',...
    'CC1_ind','alpha','rang','pair1','Npair1','sta','sta_name','tlag1',...
    'Fs','wave','vpvs','-v7.3')
