% Function to compute the mapping variable between the data and its
% backprojected image
% 
% Needs: sam_diff_tt.m, resample_tt.m
% Inputs:
%   x_no,y_no,z_no: number of nodes in x, y, and z dimensions
%   pair1,Npair1: pairs and number of station pairs (from cons_pair1.m)
%   xaxis,yaxis,zaxis: axes with coordinates along x, y, and z dimensions
%   in km
%   sta: local station coordinates in km (x,y,z)
%   Fs: sampling frequency in Hz
%   tlag1: lags corresponding to the envelope xcorr (in samples)
%   path1: path to the folder containing the travel-time tables
%   sta_name: station names to search for travel-time tables (same order as
%   in sta)
%   wave: wave type to adjust the travel-time tables
%   vpvs: needed if wave='S'
% 
% Outputs:
%   CC1_ind: sample indexes for backprojection in 3D for all station pairs
%   alpha: coeff to scale PDFs in 3D for all station pairs
% 

function [CC1_ind,alpha] = mapping_var_tt(x_no,y_no,z_no,pair1,Npair1,xaxis,...
    yaxis,zaxis,sta,Fs,tlag1,path1,sta_name,wave,vpvs)

CC1_ind = zeros(z_no,y_no,x_no,Npair1);
alpha = zeros(z_no,y_no,x_no,Npair1);

[Xg,Yg,Zg] = meshgrid(xaxis,yaxis,zaxis);
Xg = permute(Xg,[3 1 2]);
Yg = permute(Yg,[3 1 2]);
Zg = permute(Zg,[3 1 2]);

Xg = reshape(Xg,[numel(Xg) 1]);
Yg = reshape(Yg,[numel(Yg) 1]);
Zg = reshape(Zg,[numel(Zg) 1]);

dx = xaxis(2) - xaxis(1); dy = yaxis(2) - yaxis(1); dz = zaxis(2) - zaxis(1);

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
    
    % Sample difference for all grid nodes
    if strcmp(wave,'P') == 1 
        filename1 = dir([path1 '*.P.' sta_name{pair1(ipair,1),1} '.time.hdr']);
        filename1 = filename1.name(1:end-4);
        filename2 = dir([path1 '*.P.' sta_name{pair1(ipair,2),1} '.time.hdr']);
        filename2 = filename2.name(1:end-4);
        [~,ttsta1,~,~] = read_tt(path1,filename1);
        [hdrgrid,ttsta2,~,~] = read_tt(path1,filename2);
    end
    if strcmp(wave,'S') == 1 && vpvs ~= 0
        filename1 = dir([path1 '*.P.' sta_name{pair1(ipair,1),1} '.time.hdr']);
        filename1 = filename1.name(1:end-4);
        filename2 = dir([path1 '*.P.' sta_name{pair1(ipair,2),1} '.time.hdr']);
        filename2 = filename2.name(1:end-4);
        [~,ttsta1,~,~] = read_tt(path1,filename1);
        [hdrgrid,ttsta2,~,~] = read_tt(path1,filename2);
        ttsta1 = ttsta1*vpvs; ttsta2 = ttsta2*vpvs;
    end
    if strcmp(wave,'S') == 1 && vpvs == 0
        filename1 = dir([path1 '*.S.' sta_name{pair1(ipair,1),1} '.time.hdr']);
        filename1 = filename1.name(1:end-4);
        filename2 = dir([path1 '*.S.' sta_name{pair1(ipair,2),1} '.time.hdr']);
        filename2 = filename2.name(1:end-4);
        [~,ttsta1,~,~] = read_tt(path1,filename1);
        [hdrgrid,ttsta2,~,~] = read_tt(path1,filename2);
    end
    sam = sam_diff_tt( ttsta1, ttsta2, Fs );
    
    % Resample grid of travel-times if needed (and if possible)
    % Assuming dx=dy=dz
    if dx~=hdrgrid(7)
        sam = resample_tt(sam,hdrgrid,dx,max(xaxis),max(yaxis),max(zaxis));
    end

    for x_a = 1:x_no
        for y_a = 1:y_no
            for z_a = 1:z_no

                CC1_ind(z_a,y_a,x_a,ipair) = find(tlag1==sam(x_a,y_a,z_a));
                
            end
        end
    end

    clear env1 delta sta_mid_pt x_a y_a z_a
    clear s d1 d2 u sam ttsta* filename*
end

