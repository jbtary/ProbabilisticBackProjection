% Function to read the velocity model created using NonLinLoc
% 
% Input:
%   path1: path to the location of the velocity model
%   filename: filename of the velocity model without the extension
% 
% Output:
%   hdrgrid: grid paramaters (nx,ny,nz,xshift,yshift,zshift,dx,dy,dz)
%   vel: velocity grid with dims (nx,ny,nz)
% 

function [hdrgrid,vel] = read_vel(path1,filename)

% Read the hdr file of the table
fileID = fopen([path1 filename '.hdr']);
gridparams = textscan(fileID,'%d %d %d %f %f %f %f %f %f %s',1);
fclose(fileID);

% Read the buf file of the table
fileID = fopen([path1 filename '.buf']);
gridtot = fread(fileID,'float');
fclose(fileID);

% This should be the right order (checked)
vel = reshape(gridtot,[gridparams{3} gridparams{2} gridparams{1}]);
vel = permute(vel,[3 2 1]);

hdrgrid = [double(gridparams{1}) double(gridparams{2}) double(gridparams{3}) ...
    double(gridparams{4}) double(gridparams{5}) double(gridparams{6}) ...
    double(gridparams{7}) double(gridparams{8}) double(gridparams{9})];

% Velocity model often given in slowness * dc (by design in this case
% dx=dy=dz)
if strcmp(gridparams{10},'SLOW_LEN') == 1
    vel = vel/hdrgrid(7);
    vel = 1./vel;
end

% In case the velocity model is given in m/s
if max(vel(:)) > 1000
    vel = vel/1000;
end
