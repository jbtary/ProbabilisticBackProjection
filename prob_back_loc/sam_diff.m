function [ sample_diff ] = sam_diff( pair, ipair, x_a, y_a, sta, v_assume, Fs, order )
% =========================================================================
% Function: sam_diff.m
%
% Compute the sample differences between ref station and every other 
% station in the group
% =========================================================================
count = 1;

for i = pair(ipair,:)
    dist_a(count) = sqrt( (x_a - sta(i,1)).^2 + (y_a - sta(i,2)).^2 );
    count = count + 1;
end

for i = 2:order+1
    dist_diff_a(i-1) = abs( dist_a(1) - dist_a(i) );
    
    if dist_a(1) < dist_a(i)
        dist_diff_a(i-1) = -dist_diff_a(i-1);
    end
end

t_diff_a = dist_diff_a./v_assume;
sample_diff = round(t_diff_a.*Fs);