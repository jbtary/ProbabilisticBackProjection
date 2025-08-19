% Synthetic testing with different station configurations of the JC114
% dataset
clear
addpath('prob_back_loc')

% Pre-load the mapping operator and the needed variable to calculate it
% Grid parameters
grid0 = [3.0436111 -84]; % Define the 0,0 of the grid (lower-left corner, lat long)
[minx,miny] = deg2utm(grid0(1),grid0(2));
minx = minx/1000; miny = miny/1000;
rang = [0 50 0.2;0 50 0.2;0 3.2 0.2]; % Search range -- min, max, increments in km

calc = 0; % To calculate the mapping operator (=1), else you have it and load it
name_mapping = 'mapping_operator_P_whales';
name_func = 'making_distri_trans_func';

% Station configuration
stations = importdata('OBS_positions.txt',' ');
sta_name = stations.textdata;
sta_geo = [stations.data(:,2) stations.data(:,1) stations.data(:,3)]; % Lon Lat depth
[sta(:,1),sta(:,2)] = deg2utm(sta_geo(:,2),sta_geo(:,1)); sta = sta/1000;
sta(:,1) = sta(:,1) - minx;
sta(:,2) = sta(:,2) - miny;
sta(:,3) = -sta_geo(:,3); % Depth positive
Nsta = size(sta,1);

% Select the stations you want (indexes in station variable used to get
% mapping operator). Needs to be in increasing numbers.
sta_id = [1:25 29]; % [1:size(sta,1)]

% Making the mapping operator (should be done for all stations)
x_no = (rang(1,2) - rang(1,1))/rang(1,3) + 1;
y_no = (rang(2,2) - rang(2,1))/rang(2,3) + 1;
z_no = (rang(3,2) - rang(3,1))/rang(3,3) + 1;
if calc == 1
    % Calculate the mapping operator (for all stations!)
    disp('See build_mapping_operator_tt.m')
else
    % Load the mapping operator if you already have it
    load(name_mapping,'CC1_ind','alpha','pair1','Npair1')
end

% Load the transformation function from env to proba variables
%load(name_func,'bins','trans_func_all')
bins = 0.01:0.01:1; trans_func_all = make_distri_stat(bins,0.9,-9,0.001);

%% Calculate synthetics and location
clearvars -except sta z_no y_no x_no pair1 sta_id CC1_ind alpha rang ...
    bins trans_func_all sta_name
tic

% Location parameters
fact = 2; % Factor to divide PDF to calculate uncertainties
thres_prob = 1; % Threshold to reject or not a PDF from a specific pair
env_thr = 30; % Number of highest amplitude envelopes to stack

% Select the source position within rang
x_a = 10; y_a = 10; z_a = 1;

% Data and additional parameters
Fs = 100; % Sampling frequency in Hz
lag = 1900; % Max time lag included in 1st xcorr (sample) <- check with sam_a
func = 1; sig = 25; % Standard deviation of gaussian curve for env timings (in samples)
%func = 2; sig = 0.2; % Parameters for the sinc function
unc = 0.1; % If ~= 0, random time shifts to apply to timings in sec.
noise_std = 0.2; % If ~= 0, standard deviation of random noise to add to synthetics

Nsta = size(sta,1);
[pair2,Npair2] = cons_pair1(Nsta);

% Make the synthetic diff envelopes (for no transformation use trans_func_all = bins)
path1 = 'TTimes/';
rang_tt = [0 60 0.1;0 60 0.1;0 3.5 0.1]; % Grid parameters of travel-time tables (should have same 0,0 as location grid)
[env_trans,tlag1,sam_a,~] = makesyn_tt(sta,sta_name,path1,'P',rang_tt,lag,Fs,...
    x_a,y_a,z_a,func,sig,unc,noise_std,bins,trans_func_all);

max_env = max(env_trans,[],2);

CC1 = zeros(z_no,y_no,x_no,Npair2);
for ipair = 1:Npair2
    env1 = env_trans( pair2(ipair,1)==pair2(:,1) & pair2(ipair,2)==pair2(:,2), :);
    ind(ipair) = find(pair1(:,1) == sta_id(pair2(ipair,1)) & pair1(:,2) == sta_id(pair2(ipair,2)));

    indexes = squeeze(CC1_ind(:,:,:,ind(ipair)));
    CC1(:,:,:,ipair) = env1(indexes);
    %CC1(:,:,:,ipair) = CC1(:,:,:,ipair)./abs(alpha(:,:,:,ind(ipair)));

    clear indexes env1
end

prod_prob = zeros(z_no,y_no,x_no);
num_env = 0;
[~,id_sel] = maxk(max_env,env_thr);
for i = id_sel(:)' %1:Npair2
    if max_env(i) < thres_prob; continue; end
    prod_prob = prod_prob + log(CC1(:,:,:,i) + 1e-15); % Small number to avoid -Inf values (log(0))
    num_env = num_env + 1;
end
prod_prob = (1/num_env) * prod_prob;

% Get final location of event
[~,I] = max(prod_prob(:));
[iz,iy,ix] = ind2sub(size(prod_prob),I); clear I
maxp = [(ix-1)*rang(1,3) (iy-1)*rang(2,3) (iz-1)*rang(3,3)]; % Max proba position
disp(['Event located at x: ' num2str(maxp(1)) ', y: ' num2str(maxp(2)) ', z: ' num2str(maxp(3))])
clear iz iy iz

% Uncertainties (use with caution)
[uncs,uncs_coord] = calc_unc_3d(prod_prob,fact,rang);

% Calculate the RMS and residuals
%[rms_a,residuals] = calc_rms(sam_a,CC1_ind(:,:,:,ind),prod_prob,tlag1,Fs);

% Plotting
plot_env(env_trans,tlag1,pair2,id_sel(:)')
mytitle = ['Synth for source at ' num2str(x_a) ',' num2str(y_a) ',' num2str(z_a) ' km'];
plot_bkproj_prob_3d( rang, prod_prob, sta, mytitle,0)
toc
