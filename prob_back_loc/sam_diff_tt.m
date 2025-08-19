function [ sample_diff ] = sam_diff_tt( ttsta1, ttsta2, Fs )
% =========================================================================
% Function: sam_diff_tt.m
% ttsta1, ttsta2: travel-time tables for sta1 and sta2
%
% Compute the sample differences between ref station and every other 
% station in the group, using travel-time tables, for all nodes at once
% =========================================================================

% Travel-time difference at all grid nodes
if size(ttsta1) ~= size(ttsta2); disp('3D grids don t have the same size');
    sample_diff = []; return; end

diff_tt = ttsta1 - ttsta2;
sample_diff = round(diff_tt*Fs);

