% Make a function to use with back-projection location procedure
% 
% IN
%   bins: vector with values to evaluate the transfer function
%   a, b: factors to use with the exponential function
%   level: minimum value for the exponential function
% 
% OUT
%   trans_func_all: transfer function evaluated at bins values
%   curve1: curve of probability values vs envelope values
% 

function [trans_func_all,curve1] = make_distri_stat(bins,a,b,level)

curve1 = a*exp(b.*bins);
curve1(curve1<level) = level;
box = ones(1,length(bins));
box = box/sum(box);
conv_result = conv(curve1, box); % Eq. 7 numerator
trans_func_all = conv_result(1:length(bins))./curve1;

