% Build statistical distribution of xcorr values
clear
addpath(genpath('GISMO')) % Needs GISMO toolbox to load mseed data files
% Load events
path1 = 'events_picked/';
filesmat = dir([path1 '*.mat']);
ct = 1;
for ii = 1:length(filesmat)
    load([path1 filesmat(ii).name],'hdr','dath')
    dt = hdr.st(1);
    for jj = 1:length(hdr.tbeg); comp{jj} = [hdr.array{jj} hdr.comp{jj}]; end
    
    for jj = 1:length(hdr.ppicks)
        begg = round(((hdr.ppicks(jj) - hdr.tbeg(1))*86400)/dt) - (5/dt);
        endd = round(((hdr.ppicks(jj) - hdr.tbeg(1))*86400)/dt) + (15/dt);
        if begg < 0 || endd > length(dath{1}); continue; end
        
        [~,ia,~] = intersect(comp,hdr.stap(jj,:));
        datatmp = dath{ia,1};

        A2(:,ct) = datatmp(begg+1:endd);
        ct = ct+1;
    end

    clear hdr dath ia begg endd comp jj 
end
dt = 1/100;
[b,a]=butter(4,[10 45]/(1/(2*dt)),'bandpass');

% Decimate and filter
for ii = 1:size(A2,2)
    % Decimate
    A(:,ii) = decimate(A2(:,ii),500/100);
    % Filter
    A(:,ii) = filtfilt(b,a,A(:,ii));
    A(:,ii) = A(:,ii) - mean(A(:,ii));
end
clear A2 datatmp 

%% Path to continuous data
clear d_obs_*
% Stations without abnormal noises
stan = {'VA064','NG044','NG040','NG063','NG072','NG051','NG066','NG069'};
pathcont = 'JC114/DataNG';
% Get cluster info
NumTries = round(5000/size(A,2)); % Number of noise parts randomly chosen to be compared with events
[b,a]=butter(6,[10 45]/(1/(2*dt)),'bandpass');

ct = 1;
for nn = 1:length(stan)
    for day = 27:32
        fileseed = dir([pathcont '/' num2str(day) '/*' stan{nn}(end-2:end) '*hyd*.mseed']);

        % Load continuous data
        % Get random file <- Not in use
        % fid = randi(length(fileseed));
        w = waveform([pathcont '/' num2str(day) '/' fileseed(1).name],'seed');
        data = get(w,'data');
        fs = get(w,'freq');
        [p,q] = rat(dt/(1/fs));
        data = resample(data,q,p); clear p q fs w
        data = data - mean(data);
        data = filtfilt(b,a,data);
        data = data / max(abs(data));

        % Do correlation between events and random pieces of data
        kk = 1;
        for ii = 1:size(A,2)
            for jj = 1:NumTries
                % Get random piece of data
                sid = randi(length(data) - size(A,1));
                datmp1 = data(sid:sid+size(A,1)-1);
                sid = randi(length(data) - size(A,1));
                datmp2 = data(sid:sid+size(A,1)-1);
                c_noise(kk,:) = abs(hilbert(xcorr(datmp1,datmp2,'coeff')));
                c_sig(kk,:) = abs(hilbert(xcorr(A(:,ii),datmp1,'coeff')));
                kk = kk+1;
                clear sid datmp*
            end
        end
        
        binedges = 0:0.01:1;
        d_obs_noise(:,ct) = histcounts(c_noise(:),binedges, 'Normalization', 'probability');
        d_obs_sig(:,ct) = histcounts(c_sig(:),binedges, 'Normalization', 'probability');
        ct = ct+1;
        clear c_* data
    end
end
clear fileseed ii jj nn 
%save making_distri_trans_func

%%
bins = 0.005:0.01:1;
d_obs_noise_m = mean(d_obs_noise,2);
d_obs_sig_m = mean(d_obs_sig,2);
figure; 
subplot(131); plot(bins,d_obs_noise)
subplot(132); plot(bins,d_obs_sig)
subplot(133); plot(bins,d_obs_noise_m,'b-*',bins,d_obs_sig_m,'k-*')

% Fit distribution (both exp and gauss work really well)
%ft = fittype('a*exp(b*x)','independent','x','dependent','y');
%ft = fittype('a*exp(-((x-b)/c)^2)','independent','x','dependent','y');
%opts = fitoptions('Method','NonlinearLeastSquares');
%opts.StartPoint = [0.1 0.1 0.1];
%opts.Lower = [-Inf -Inf 0];
%opts.Upper = [Inf Inf Inf];
%fit_obj1 = fit(bins',d_obs_noise_m,ft,opts);

% distri_a = fit_rat1_a*exp(-fit_rat1_a*bins);

level = 0.0001;
d_obs_noise_m(d_obs_noise_m<level) = level;
d_obs_sig_m(d_obs_sig_m<level) = level;

box = ones(1,length(bins));
box = box/sum(box);
conv_result = conv(d_obs_noise_m, box); % Eq. 7 numerator
trans_func_all = conv_result(1:length(bins))./d_obs_sig_m';
figure; subplot(131); plot(bins,d_obs_noise_m,'b-*',bins,d_obs_sig_m,'k-*')
subplot(132); plot(conv_result)
subplot(133); plot(bins,trans_func_all)

%save making_distri_trans_func