function [ pair1, Npair1 ] = cons_pair1( Nsta )
% =========================================================================
% Function: cons_pair1.m
%
% Compute all combinations of station pairs (1st order)
%
% Number of combination: nCk     n = total no. of station
%                                k = no. of station in a group
% =========================================================================
count = 1;
for i = 1:Nsta
    for j = 1:Nsta
        if j > i
            pair1(count,:) = [i j];
            count = count + 1;
        end
    end
end

Npair1 = length(pair1(:,1));