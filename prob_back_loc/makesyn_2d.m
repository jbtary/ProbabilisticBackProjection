% Make synthetic data for tremor location, for a 2D constant velocity
% structure.
% 
% Needs: cons_pair1.m
% Input:
%   sta: station positions in x,y (a priori in km)
%   vel: constant velocity (a priori in km/s)
%   lag: max lag to consider in samples (must be larger than max lag between stations)
%   Fs: sampling frequency in Hz
%   x_a,y_a: position of the fake source (same unit as station coordinates)
%   func: number to choose env function shape (1: gaussian, 2: sinc)
%   sig: parameters needed for each function "func" (sigma for gaussian, main lobe width 
%   for sinc)
%   unc: random time shifts to apply to the timings in sec. (=0 no time shift, else >0)
%   nstd: standard deviation for some random noise 
% 
% Output:
%   env_trans: synthetic differential envelopes
%   tlag1: lag axis
%   sam_a: differential timings in samples for each station pair
%   alpha: amplitude scaling parameter

function [env_trans,tlag1,sam_a,alpha] = makesyn_2d(sta,vel,lag,Fs,x_a,y_a,func,sig,unc,nstd)

Nsta = size(sta,1);
[pair1,Npair1] = cons_pair1(Nsta);
tlag1 = -lag:lag;

for ii = 1:size(pair1,1)
    % 2D euclidean distance stations - source
    dist_a(1) = norm([x_a y_a] - sta(pair1(ii,1),:));
    dist_a(2) = norm([x_a y_a] - sta(pair1(ii,2),:));
    dist_diff_a = abs( dist_a(1) - dist_a(2) );
    if dist_a(1) < dist_a(2); dist_diff_a = -dist_diff_a; end
    % Fake xcorr lag data: 1st step
    sam_a(ii) = round(dist_diff_a/vel*Fs); % Corresponding sample timing
    
    % Add a random time shift to represent uncertainties (time and/or velocity model)
    if unc ~= 0; sam_a(ii) = sam_a(ii) + round(((rand(1)*unc)-unc/2)*Fs); end

    % Calculation of alpha
    delta = norm(sta(pair1(ii,2),:) - sta(pair1(ii,1),:)); % distance between stations
    sta_mid_pt = [(sta(pair1(ii,1),1)+sta(pair1(ii,2),1))/2 ...
        (sta(pair1(ii,1),2)+sta(pair1(ii,2),2))/2]; % mid-point coordinates
    s = norm([x_a y_a] - sta_mid_pt); % distance source - mid point
    u = dist_a(1) - dist_a(2); % diff of distances
    alpha(ii) = sqrt((s^2+delta^2/4-u^2/2) / (delta^2-u^2));
end

if any(isinf(alpha))
    disp('Some alpha values were Inf') % When u and s are the same
    alpha(isinf(alpha)) = 1; 
end

if abs(max(sam_a)) > abs(max(tlag1)); disp('Max lag too small'); env_trans = []; return; end

env1 = zeros(Npair1,length(tlag1));

for ii = 1:size(pair1,1)
    env1(ii,tlag1==sam_a(ii)) = 1; % Add fake sample diff to each pair: 2nd step
end

% Gaussian curve
gaus = @(x,mu,sig,amp)amp*exp(-(((x-mu).^2)/(2*sig.^2)));
if func == 1 
    x = -sig*5:sig*5; mu = 0; amp = 1;
    y = gaus(x,mu,sig,amp);
end
% Sinc function
if func == 2
    x = linspace(-10/sig,10/sig,round(100/sig));
    y = sinc(sig*x); % sig>1 to have thinner main lobe, 0<sig<1 to have wider ones
    % Multiply by the inverse of a Hann window to boost amplitudes of side-lobes
    win = 1.5*(1-gaus(x,0,1/sig,1)) + 0.5;
    y = y.*win; 
    win = tukeywin(length(y),0.05);
    y = abs(y.*win');
end 

y = y/max(y);
env_trans = zeros(Npair1,length(tlag1));
for ii = 1:size(pair1,1)
    env_trans(ii,:) = conv(env1(ii,:),y,'same'); % conv with gaussian: 3rd step
    % Add some random noise
    if nstd > 0
        noise = nstd*abs(randn(size(tlag1))); noise = medfilt1(noise,20);
        env_trans(ii,:) = env_trans(ii,:) + noise;
        env_trans(ii,:) = env_trans(ii,:)/max(env_trans(ii,:));
    end
    env_trans(ii,:) = env_trans(ii,:)*alpha(ii);
end
