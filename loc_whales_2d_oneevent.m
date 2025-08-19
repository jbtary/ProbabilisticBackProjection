% Location of events using the method and codes of Ka Lok Li for tremors
% Version to locate 1 event at a time in 2d
clear; close all
addpath('prob_back_loc')

% Get the whale files
allsta = 'OBS_positions.txt'; % File with all station names and locations

% Load precalculated mapping operator (all other variables assumed to be
% the same)
name_mapping = 'mapping_operator_whales_2d';
load(name_mapping,'CC1_ind','alpha','max_samlag','pair1','x_no','y_no')
% Load or make the function to transform CC envelope values to
% probabilities
%name_trans_func = 'making_distri_trans_func';
%load(name_trans_func,'bins','trans_func_all'); bin_val = bins; clear bins
bin_val = 0.01:0.01:1; trans_func_all = make_distri_stat(bin_val,0.9,-9,0.001);

%% Other parameters to set
pathtoev = 'Events_mat/';
files = dir([pathtoev '*.mat']); 
filenum = 1;
grid0 = [3.0436111 -84.0]; % 0,0 of grids (lower-left corner, lat long)
% lag must be > all max_samlag and < length(sig)
lag = 1900; % Max time lag included in 1st xcorr (sample)
fact = 2; % Factor to divide PDF to calculate uncertainties
thres_prob = 1; % Threshold to reject or not a PDF from a specific pair
env_thr = 30; % Number of highest amplitude envelopes to stack
svpath = 'Save/';

% Study area and stations
rang = [0 50 0.2;0 50 0.2];         % Search range -- min, max, increment in km
% Filtering (Butterworth)
Fs = 100;                           % Sampling frequency for re-sampling in Hz
freq1 = 10;                         % Lower corner freq (Hz)
freq2 = 45;                         % Upper corner freq (Hz)
filt_order = 4;                     % Butterworth filter order
% One-bit normalization
one_bit_opt = 0;                    % Doing one-bit normalization? 1=yes, 0=no

if ~exist(svpath,'dir'); mkdir(svpath); end

% Cleaning
clearvars -except CC1_ind alpha max_samlag pair1 x_no y_no grid0 ...
    pathtoev lag fact jj files allsta svpath ...
    Loct_sv name_* rang Fs freq1 freq2 filt_order one_bit_opt ...
    trans_func_all bin_val filenum thres_prob env_thr

tic
disp(['File ' files(filenum).name])
load([pathtoev files(filenum).name],'hdr')

% Select stations (all: sta_sel = {};)
ct = 1;
for ii = 1:size(hdr.comp,2)
    if strcmp('NG',hdr.array{ii}); sta_sel{ct} = [hdr.array{ii} hdr.comp{ii}]; ct = ct+1; end
    if strcmp('VA364',[hdr.array{ii} hdr.comp{ii}]); sta_sel{ct} = 'VA364'; ct = ct+1; end
end

[sta_name,sta,ind,data,hdr2] = sta_extract_h([pathtoev files(filenum).name],allsta,grid0,sta_sel);

% Define a timing for the event inside of file in seconds
win = [0 30];

% Check if all data streams are at least "win" long
len_d = cellfun(@length,data.h);
ind_d = (len_d*hdr2.st(1) >= win(2) - win(1));
if sum(ind_d == 0) > 0
    data.h = data.h(ind_d == 1);
    hdr2.st = hdr2.st(ind_d == 1); hdr2.tbeg = hdr2.tbeg(ind_d == 1);
    sta_name = sta_name(ind_d == 1);
    sta = sta(ind_d == 1);
    ind = ind(ind_d == 1);
end
clear len_d ind_d

% -------------------------------------------------------------------------
% Read data and processing
% -------------------------------------------------------------------------
disp('Data processing')

% Read data files
Nsta = size(sta,1); % Number of stations
for ii = 1:Nsta
    sig1 = data.h{ii,1};

    % Get the date of seismogram (same for all data)
    startdate = hdr2.tbeg(ii);
    % Get the sampling rate (Hz) (same for all data)
    if length(hdr2.st) > 1
        fs1 = round(1/hdr2.st(ii));
    else
        fs1 = round(1/hdr2.st);
    end

    sig1 = decimate(sig1,fs1/Fs);

    % ---------------------------------------------------------------------
    % Extract seismogram corresponding to current time window
    % ---------------------------------------------------------------------
    sig2(ii,:) = sig1(round(win(1)*Fs+1) : round(win(2)*Fs));
    t_axis_seg(ii,:) = win(1)+(1/Fs):(1/Fs):win(2);

    sig2(ii,:) = sig2(ii,:) - mean(sig2(ii,:));

    clear sig1 fs1
end

% Filter the signals
for ii = 1:Nsta
    sig3(ii,:) = butterworthFilt(sig2(ii,:)','bandpass',filt_order,Fs,freq1,freq2);
    sig3(ii,:) = sig3(ii,:) * 42.6e-6; % Total sensitivity in Pa/count
end
clear sig2

% One-bit normalization
if one_bit_opt == 1
    for ii = 1:Nsta
        sig(ii,:) = sign(sig3(ii,:));
    end
elseif one_bit_opt == 0
    sig = sig3;
end
clear sig3

% ------------------------------------------------------------------------
% Cross-correlations and its envelopes
% -------------------------------------------------------------------------
disp('CC and envelopes')
[ env, tlag1, pairev ] = FCC_all_3d( Nsta, sig, lag );

% Select the pairs needed depending on the station in current event
for ii = 1:size(pairev,1)
    % Pair indexes in current event data followed by corresponding indexes
    % in mapping operator
    iv = find(pair1(:,1) == ind(pairev(ii,1)) & pair1(:,2) == ind(pairev(ii,2)));
    pair2(ii,:) = [pairev(ii,:) iv];
    max_samlag2(ii) = max_samlag(iv);
end
Npair2 = size(pair2,1); clear iv

% Smoothing of the envelopes using a low pass filter
env = smooth_env(env,Fs,5,2);

% -------------------------------------------------------------------------
% Doing the transformation from covariogram envelopes to probabilities
% -------------------------------------------------------------------------
disp('Transforming covariograms to probabilities')

env_trans = zeros(size(env));
for i = 1:size(env,1)
    if mod(i,100) == 0; disp(['Pair no ' num2str(i) '/' num2str(size(env,1))]); end
    for j = 1:size(env,2)
        ind2 = find(bin_val > env(i,j),1);

        if ind2 >= 2
            ind1 = ind2 - 1;
            env_trans(i,j) = trans_func_all(1,ind1) + ...
                (env(i,j) - bin_val(ind1))*(trans_func_all(1,ind2)-trans_func_all(1,ind1))...
                /(bin_val(ind2)-bin_val(ind1));
        else
            env_trans(i,j) = trans_func_all(1,ind2);
        end

    end
end
clear ind2 ind1

for ii = 1:size(env,1)
    ind_s = find(tlag1 == -max_samlag2(ii));
    ind_e = find(tlag1 == max_samlag2(ii));
    max_env(ii) = max(env_trans(ii,ind_s:ind_e));
end
clear ind_*

% -------------------------------------------------------------------------
% Compute the back projected map for each pair
% -------------------------------------------------------------------------
CC1 = zeros(y_no,x_no,Npair2);
for ipair = 1:Npair2
    env1 = env_trans( pair2(ipair,1)==pair2(:,1) & pair2(ipair,2)==pair2(:,2), :);

    indexes = squeeze(double(CC1_ind(:,:,pair2(ipair,3))));
    CC1(:,:,ipair) = env1(indexes);
    % Absolute value of alpha because it happens that alpha is complex
    % (very very rarely)
    %CC1(:,:,ipair) = CC1(:,:,ipair)./abs(alpha(:,:,pair2(ipair,3)));

    clear indexes env1
end

% -------------------------------------------------------------------------
% Calculate the product of probabilities from all station pairs
% -------------------------------------------------------------------------
prod_prob = zeros(y_no,x_no);
num_env = 0;
[~,id_sel] = maxk(max_env,env_thr);
for i = id_sel %1:Npair2
    if max_env(i) < thres_prob; continue; end
    prod_prob = prod_prob + log(CC1(:,:,i));
    num_env = num_env + 1; % Quality measure of PDF as number of pairs contributing to solution
end

prod_prob = (1/num_env) * prod_prob;

% Uncertainties (use with caution)
[uncs,uncs_coord] = calc_unc_2d(prod_prob,fact,rang);

% Get max of probability volume
[~,I] = max(prod_prob(:));
[iy,ix] = ind2sub(size(prod_prob),I); clear I
maxp = [(ix-1)*rang(1,3) (iy-1)*rang(2,3)]; % Max proba position
[minx,miny,utmzone] = deg2utm(grid0(1),grid0(2));
utm = [(maxp(1)*1000) + minx (maxp(2)*1000) + miny];
[lat,lon] = utm2deg(utm(1),utm(2),utmzone);

% YYYY MM DD HH MM SS.SSS LAT LON DEPTH MAXP NPH GAP DIST ERRH ERRZ
tim_st = datestr(hdr2.tbeg(1),'dd/mm/yyyy HH:MM:SS.FFF');
Loct = [str2double(tim_st(7:10)) str2double(tim_st(4:5)) str2double(tim_st(1:2)) ...
    str2double(tim_st(12:13)) str2double(tim_st(15:16)) str2double(tim_st(18:end)) ...
    lat lon 0 max(prod_prob(:)) num_env 0 0 max(uncs) 0];

% ---------------------------------------------------------------------
% Plotting
% ---------------------------------------------------------------------
plot_bkproj_prob_2d( rang, prod_prob, sta, files(filenum).name, 0)
plot_ell(maxp(1),maxp(2),uncs(1)/2,uncs(2)/2)
filename = ['Locs_pdf_' files(filenum).name(1:end-4)];
print('-dpng','-r150',[svpath filename])

save([svpath filename],'prod_prob','rang','sta','sta_name','Loct',...
    'env','env_trans','freq1','freq2','Fs','tlag1','pair2','max_env',...
    'max_samlag2') % Quite heavy...

toc

% Display location results
disp('Location results:')
disp('YYYY MM DD HH MM SS.SSS LAT LON DEPTH MAXP NPH GAP DIST ERRH ERRZ')
str = sprintf('%i %i %i %i %i %2.3f %2.6f %2.6f %i %2.2f %i %i %i %2.3f %i\n',Loct');
disp(str)
