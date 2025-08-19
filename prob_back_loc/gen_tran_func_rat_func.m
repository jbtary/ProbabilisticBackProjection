function [ trans_func_all, bin_val ] = gen_tran_func_rat_func( env_norm, tlag1, max_samlag, bin_resol )
% =========================================================================
% Function: gen_tran_func_rat_func.m
%
% Generate the function used to transform the time series of covariogram
% envelope to the time series of probabilities.
%
% In this version:
% - the distribution of noise is estimated from the empirical distribution
%   of time series (excluding amplitudes that are more than 75 % of the
%   maximum)
% - the distribution of time series (signal + noise) is estimated from the
%   whole empirical distribution of time series
% =========================================================================

% Construct a vector of bin values up to the global max. of all station
% pairs
bin_edges = (0:bin_resol:ceil(max(max(env_norm)))+bin_resol);
for i = 1:length(bin_edges)-1
    bin_val(i) = (bin_edges(i)+bin_edges(i+1))/2;
end

trans_func_all = zeros(length(env_norm(:,1)),length(bin_val));

% Loop through all station pairs
for i = 1:length(env_norm(:,1))
    clear theo_noise box theo_t_series
    
    ind_s = find(tlag1 == -max_samlag(i));
    ind_e = find(tlag1 == max_samlag(i));
    
    % The empirical distribution of covariogram envelope
    [cc,~] = histcounts(env_norm(i,ind_s:ind_e),bin_edges);
    
    % Find the normalization constant
    ss = 0;
    for j = 1:length(bin_val)
        ss = ss + cc(j)*bin_resol;
    end
    
    % Normalize distribution so that it integrates to 1
    cc = cc./ss;
    
    % Find the amplitude at 75 % of the maximum of covariogram envelope
    max_noise_amp = 0.75*max(env_norm(i,ind_s:ind_e));
    ind = find(bin_val > max_noise_amp,1);
    
    % Find the distribution of noise
    ft = fittype('(x/a)/(1+(x/a)^b)','independent','x','dependent','y');
    opts = fitoptions('Method','NonlinearLeastSquares');
    opts.StartPoint = [0.1 1];
    opts.Lower = [0.1 0.1];
    opts.Upper = [10 10];
    fit_obj1 = fit(bin_val(1:ind)',cc(1:ind)',ft,opts);
    fit_rat1_a = fit_obj1.a;
    fit_rat1_b = fit_obj1.b;
    
    % Construct the theoretical exponential distribution for noise: Eq. 8
    theo_noise = (bin_val./fit_rat1_a) ./ (1 + (bin_val./fit_rat1_a).^fit_rat1_b);
    ss = 0;
    for j = 1:length(bin_val)
        ss = ss + theo_noise(j)*bin_resol;
    end
    theo_noise = theo_noise./ss;
    
    % Find the distribution of time series (signal + noise)
    fit_obj2 = fit(bin_val',cc',ft,opts);
    fit_rat2_a = fit_obj2.a;
    fit_rat2_b = fit_obj2.b;
    
    % Construct the theoretical exponential distribution for time series
    % (signal + noise): Eq. 8
    theo_t_series = (bin_val./fit_rat2_a) ./ (1 + (bin_val./fit_rat2_a).^fit_rat2_b);
    ss = 0;
    for j = 1:length(bin_val)
        ss = ss + theo_t_series(j)*bin_resol;
    end
    theo_t_series = theo_t_series./ss;
    
    % Boxcar distribution (distribution for signals only)
    box = ones(1,length(bin_val));

    ss = 0;
    for j = 1:length(bin_val)
        ss = ss + box(j)*bin_resol;
    end
    box = box./ss;

    % Convolution of distribution of noise and signals
    conv_result = conv(theo_noise, box); % Eq. 7 numerator

    % The required transformation
    trans_func_all(i,:) = conv_result(1:length(bin_val))./theo_t_series; % Eq. 7
end
