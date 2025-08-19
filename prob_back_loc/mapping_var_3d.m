% Function to compute the mapping variable between the data and its
% backprojected image
% 
% Needs: sam_diff_3d.m
% Inputs:
%   x_no,y_no,z_no: number of nodes in x, y, and z dimensions
%   pair1,Npair1: pairs and number of station pairs (from cons_pair1.m)
%   xaxis,yaxis,zaxis: axes with coordinates along x, y, and z dimensions
%   in km
%   sta: local station coordinates in km (x,y,z)
%   v_assume: constant velocity in km/s
%   Fs: sampling frequency in Hz
%   tlag1: lags corresponding to the envelope xcorr (in samples)
% 
% Outputs:
%   CC1_ind: sample indexes for backprojection in 3D for all station pairs
%   alpha: coeff to scale PDFs in 3D for all station pairs
% 

function [CC1_ind,alpha] = mapping_var_3d(x_no,y_no,z_no,pair1,Npair1,xaxis,...
    yaxis,zaxis,sta,v_assume,Fs,tlag1)

CC1_ind = zeros(z_no,y_no,x_no,Npair1);
alpha = zeros(z_no,y_no,x_no,Npair1);

[Xg,Yg,Zg] = meshgrid(xaxis,yaxis,zaxis);
Xg = permute(Xg,[3 1 2]);
Yg = permute(Yg,[3 1 2]);
Zg = permute(Zg,[3 1 2]);

Xg = reshape(Xg,[numel(Xg) 1]);
Yg = reshape(Yg,[numel(Yg) 1]);
Zg = reshape(Zg,[numel(Zg) 1]);

for ipair = 1:Npair1
    fprintf('Processing pair %d/%d, %d   %d\n',ipair,Npair1,pair1(ipair,1),pair1(ipair,2));

    % Distance between two stations
    delta = norm(sta(pair1(ipair,2),:) - sta(pair1(ipair,1),:));

    % Coordinates of mid point between two stations
    sta_mid_pt = [(sta(pair1(ipair,1),1)+sta(pair1(ipair,2),1))/2 ...
        (sta(pair1(ipair,1),2)+sta(pair1(ipair,2),2))/2 ...
        (sta(pair1(ipair,1),3)+sta(pair1(ipair,2),3))/2];
    
    % Distance from current grid point to mid point
    s = vecnorm([Xg Yg Zg]' - sta_mid_pt');
    s = reshape(s,[z_no y_no x_no]);
    % Distance from station 1 to current grid point
    d1 = vecnorm([Xg Yg Zg]' - sta(pair1(ipair,1),:)');
    d1 = reshape(d1,[z_no y_no x_no]);
    % Distance from station 2 to current grid point
    d2 = vecnorm([Xg Yg Zg]' - sta(pair1(ipair,2),:)');
    d2 = reshape(d2,[z_no y_no x_no]);
    % Differential distance between two stations
    u = d1 - d2;
    % Factor used to compensate for hyperbolic effects: Eq. 9
    alpha(:,:,:,ipair) = sqrt( (s.^2+delta^2/4-u.^2/2) ./ (delta^2-u.^2) );

    xcount = 1;
    for x_a = xaxis
        ycount = 1;
        for y_a = yaxis
            zcount = 1;
            for z_a = zaxis

                sam = sam_diff_3d( pair1, ipair, x_a, y_a, z_a, sta, v_assume, Fs, 1 );
                CC1_ind(zcount,ycount,xcount,ipair) = find(tlag1==sam(1,1));
                
                zcount = zcount + 1;
            end
            ycount = ycount + 1;
        end
        xcount = xcount + 1;
    end

    clear env1 delta sta_mid_pt xcount ycount zcount x_a y_a z_a
    clear s d1 d2 u sam
end

