% Function to calculate the RMS from the sample position of max of
% envelope, the max of the proba position, and the mapping operator, for
% all station pairs used for the location
% 
% Input:
%   sam_a: position of max of envelope in samples, for all station pairs
%   used to do the location
%   CC1_ind: mapping operator to do the location (only for the station
%   pairs included in "sam_a")
%   prod_prob: final PDF of the location
%   tlag1: time axis in samples to search for the right lag time
%   corresponding to a diff time in CC1_ind
%   Fs: sampling frequency in Hz
% 
% Output:
%   rms_a: rms in sec.
%   res_s: residuals for each station pair in sec.

function [rms_a,res_s] = calc_rms(sam_a,CC1_ind,prod_prob,tlag1,Fs)

% Get max of probability cube
[~,I] = max(prod_prob(:));
[iz,iy,ix] = ind2sub(size(prod_prob),I);

% Get the lag time corresponding to max for each station pair
for ii = 1:length(sam_a) % Loop on the pairs of stations
    sam_m = CC1_ind(iz,iy,ix,ii);
    tlag_m = tlag1(sam_m);
    res_s(ii) = (sam_a(ii) - tlag_m)/Fs;
end

% Calculate rms
rms_a = rms(res_s);

