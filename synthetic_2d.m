% Synthetic testing with different station configurations of the JC114
% dataset
clear
addpath('prob_back_loc')

% Grid parameters
grid0 = [3.0436111 -84]; % Define the 0,0 of the grid (lower-left corner, lat long)
[minx,miny] = deg2utm(grid0(1),grid0(2));
minx = minx/1000; miny = miny/1000;
rang = [0 50 0.2;0 50 0.2]; % Search range -- min, max, increments in km

calc = 0; % To calculate the mapping operator (=1), else you have it and load it
name_mapping = 'mapping_operator_whales_2d';

% Station configuration
stations = importdata('OBS_positions.txt',' ');
sta_name = stations.textdata;
sta_geo = [stations.data(:,2) stations.data(:,1) stations.data(:,3)]; % Lon Lat depth
[sta(:,1),sta(:,2)] = deg2utm(sta_geo(:,2),sta_geo(:,1)); sta = sta/1000;
sta(:,1) = sta(:,1) - minx;
sta(:,2) = sta(:,2) - miny;
Nsta = size(sta,1);

% Select the stations you want (indexes in station variable used to get
% mapping operator). Needs to be in increasing numbers.
sta_id = [1:25 29];

% Making the mapping operator (should be done for all stations)
x_no = (rang(1,2) - rang(1,1))/rang(1,3) + 1;
y_no = (rang(2,2) - rang(2,1))/rang(2,3) + 1;
if calc == 1
    % Calculate the mapping operator (for all stations!)
    Fs = 100;
    vel = 1.5;
    [pair1,Npair1] = cons_pair1(Nsta);
    tlag1 = -lag:lag;
    [CC1_ind,alpha] = mapping_var_2d(x_no,y_no,pair1,Npair1,rang(1,1):rang(1,3):rang(1,2),...
        rang(2,1):rang(2,3):rang(2,2),sta,vel,Fs,tlag1);

    % Save the mapping operator for future uses if needed
    save(name_mapping,'x_no','y_no','CC1_ind','alpha',...
        'rang','pair1','Npair1','sta','sta_name','tlag1','Fs','vel')
else
    % Load the mapping operator if you already have it
    load(name_mapping,'CC1_ind','alpha','pair1','Npair1')
end

%% Calculate synthetics and location
clearvars -except sta y_no x_no pair1 sta_id CC1_ind alpha rang ...
    bins trans_func_all 
tic

% Location parameters
fact = 2; % Factor to divide PDF to calculate uncertainties
%thres_prob = 0; % Threshold to reject or not a PDF from a specific pair
%env_thr = 325; % Number of highest amplitude envelopes to stack
lim_uncs = 1e20; % Limit to calculate the uncertainties

% Select the source position within rang
x_a = 5; y_a = 30;

% Data and additional parameters
Fs = 100; % Sampling frequency in Hz
lag = 1900; % Max time lag included in 1st xcorr (sample)
vel = 1.5; % Back projection velocity (km/s)
func = 1; sig = 20; % Standard deviation of gaussian curve for env timings (in samples)
%func = 2; sig = 0.2; % Parameters for the sinc function
unc = 0; % If ~= 0, random time shifts to apply to timings in sec.
noise_std = 0; % If ~= 0, standard deviation of random noise to add to synthetics

Nsta = size(sta,1);
[pair2,Npair2] = cons_pair1(Nsta);

% Make the synthetic diff envelopes
[env_trans,tlag1,sam_a,~] = makesyn_2d(sta,vel,lag,Fs,x_a,y_a,func,sig,unc,noise_std);

CC1 = zeros(y_no,x_no,Npair2);
for ipair = 1:Npair2
    env1 = env_trans( pair2(ipair,1)==pair2(:,1) & pair2(ipair,2)==pair2(:,2), :);
    ind(ipair) = find(pair1(:,1) == sta_id(pair2(ipair,1)) & pair1(:,2) == sta_id(pair2(ipair,2)));
    
    indexes = squeeze(CC1_ind(:,:,ind(ipair)));
    CC1(:,:,ipair) = env1(indexes);
    CC1(:,:,ipair) = CC1(:,:,ipair)./alpha(:,:,ind(ipair));
    
    clear indexes env1
end

prod_prob = ones(y_no,x_no);
num_env = 0;
%[~,id_sel] = maxk(max_env,env_thr);
for i = 1:Npair2 %id_sel(:)'
    %if max_env(i) < thres_prob; continue; end
    prod_prob = prod_prob + log(CC1(:,:,i) + 1e-15); % Small number to avoid -Inf values (log(0))
    num_env = num_env + 1;
end
prod_prob = (1/num_env) * prod_prob;

% Get calculated location of event
[~,I] = max(prod_prob(:));
[iy,ix] = ind2sub(size(prod_prob),I); clear I
maxp = [(ix-1)*rang(1,3) (iy-1)*rang(2,3)]; % Max proba position
disp(['Event located at x: ' num2str(maxp(1)) ', y: ' num2str(maxp(2))])

% Calculate the RMS and residuals
for ii = 1:length(sam_a) % Loop on the pairs of stations
    sam_m = CC1_ind(iy,ix,ind(ii));
    tlag_m = tlag1(sam_m);
    residuals(ii) = (sam_a(ii) - tlag_m)/Fs;
end
clear sam_m tlag_m
rms_a = rms(residuals);

% Uncertainties (use with caution)
[uncs,uncs_coord] = calc_unc_2d(prod_prob,fact,rang);

% Plotting
plot_env(env_trans,tlag1,pair2,1:10)
mytitle = ['Synth for source at ' num2str(x_a) ',' num2str(y_a) ' km'];
plot_bkproj_prob_2d( rang, prod_prob, sta, mytitle, 0 )
toc
