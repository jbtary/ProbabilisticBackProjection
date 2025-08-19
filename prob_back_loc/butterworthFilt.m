function [sig_out] = butterworthFilt(sig,ftype,n,fs,f1,f2)
% -----------------------------------------------------------------------
% function: butterworthFilt.m
% Written by Ka Lok Li, Uppsala University
% Last modified: 2015/03/26
%
% Filter a time series using the Butterworth filter
%
% Input:
% sig   = input time series 
% ftype = filter type ( 'low' | 'bandpass' | 'high' | 'stop' )
% n     = filter order
% fs    = sampling frequency (Hz)
% f1    = cutoff frequency 1 (Hz)
% f2    = cutoff frequency 2 (Hz) (optional)
%
% Output:
% sig_out = filtered time series
% -----------------------------------------------------------------------
% Nyquist frequency
fnq = fs/2;

% Cutoff frequency
if nargin == 6
    Wn = [f1/fnq f2/fnq];
elseif nargin == 5
    Wn = f1/fnq;
end

% Design a Butterworth filter 
[z,p,k] = butter(n,Wn,ftype); 
[sos,g] = zp2sos(z,p,k);
sig_out = filtfilt(sos,g,sig);

end
