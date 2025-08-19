% Smooth the envelopes using a low-pass zero-phase filter
% 
% IN
%   env: envelopes to be filtered (Nsta x Nsamples)
%   fs: sampling frequency in Hz
%   fmax: maximum frequency of low pass filter in Hz
%   fo: filter order (ex: 2)
% OUT
%   envf: filtered envelopes (Nsta x Nsamples)

function envf = smooth_env(env,fs,fmax,fo)

max_env = max(env,[],2);
[b,a] = butter(fo,fmax/(fs/2),'low'); % Filter parameters
envf = filtfilt(b,a,env'); envf = envf';
max_env2 = max(envf,[],2);
rat = max_env./max_env2;
rat_mat = repmat(rat,1,size(env,2));
envf = envf.*rat_mat;

