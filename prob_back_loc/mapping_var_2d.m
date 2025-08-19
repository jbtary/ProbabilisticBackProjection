% Function to compute the mapping variable between the data and its
% backprojected image
% 
% Needs: sam_diff.m
% Inputs:
%   x_no,y_no: number of nodes in x, and y dimensions
%   pair1,Npair1: pairs and number of station pairs (from cons_pair1.m)
%   xaxis,yaxis: axes with coordinates along x, and y dimensions in km
%   sta: local station coordinates in km (x,y)
%   v_assume: constant velocity in km/s
%   Fs: sampling frequency in Hz
%   tlag1: lags corresponding to the envelope xcorr (in samples)
% 
% Outputs:
%   CC1_ind: sample indexes for backprojection in 2D for all station pairs
%   alpha: coeff to scale PDFs in 2D for all station pairs
% 

function [CC1_ind,alpha] = mapping_var_2d(x_no,y_no,pair1,Npair1,xaxis,...
    yaxis,sta,v_assume,Fs,tlag1)

CC1_ind = zeros(y_no,x_no,Npair1);
alpha = zeros(y_no,x_no,Npair1);

[Xg,Yg] = meshgrid(xaxis,yaxis); % Matrices y_no*x_no

Xg = reshape(Xg,[numel(Xg) 1]);
Yg = reshape(Yg,[numel(Yg) 1]);

for ipair = 1:Npair1
    fprintf('Processing pair %d/%d, %d   %d\n',ipair,Npair1,pair1(ipair,1),pair1(ipair,2));

    % Distance between two stations
    delta = norm(sta(pair1(ipair,2),:) - sta(pair1(ipair,1),:));

    % Coordinates of mid point between two stations
    sta_mid_pt = [(sta(pair1(ipair,1),1)+sta(pair1(ipair,2),1))/2 ...
        (sta(pair1(ipair,1),2)+sta(pair1(ipair,2),2))/2];
    
    % Distance from current grid point to mid point
    s = vecnorm([Xg Yg]' - sta_mid_pt');
    s = reshape(s,[y_no x_no]);
    % Distance from station 1 to current grid point
    d1 = vecnorm([Xg Yg]' - sta(pair1(ipair,1),:)');
    d1 = reshape(d1,[y_no x_no]);
    % Distance from station 2 to current grid point
    d2 = vecnorm([Xg Yg]' - sta(pair1(ipair,2),:)');
    d2 = reshape(d2,[y_no x_no]);
    % Differential distance between two stations
    u = d1 - d2;
    % Factor used to compensate for hyperbolic effects: Eq. 9
    alpha(:,:,ipair) = sqrt( (s.^2+(delta^2/4)-(u.^2/2)) ./ (delta^2-u.^2) );

    xcount = 1;
    for x_a = xaxis
        ycount = 1;
        for y_a = yaxis

            sam = sam_diff( pair1, ipair, x_a, y_a, sta, v_assume, Fs, 1 );
            CC1_ind(ycount,xcount,ipair) = find(tlag1==sam(1,1));
            
            ycount = ycount + 1;
        end
        xcount = xcount + 1;
    end

    clear env1 delta sta_mid_pt xcount ycount x_a y_a
    clear s d1 d2 u sam
end

