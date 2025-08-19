% Function to read the travel-time tables created using NonLinLoc
% 
% Input:
%   path1: path to the location of the travel-time tables
%   filename: filename of the travel-time table without the extension
% 
% Output:
%   hdrgrid: grid paramaters (nx,ny,nz,xshift,yshift,zshift,dx,dy,dz)
%   tt: travel time grid with dims (nx,ny,nz)
%   staparams: station parameters (name,x,y,z)
%   grid_ori: grid origin (lat,long of lower left corner)
% 

function [hdrgrid,tt,staparams,grid_ori] = read_tt(path1,filename)

% Read the hdr file of the table
fileID = fopen([path1 filename '.hdr']);
gridparams = textscan(fileID,'%d %d %d %f %f %f %f %f %f %s %s',1);
staparams = textscan(fileID,'%s %f %f %f',1);
gridori = textscan(fileID,'%s %s %s %f %s %f %s %f',1);
fclose(fileID);

% Read the buf file of the table
fileID = fopen([path1 filename '.buf']);
gridtot = fread(fileID,'float');
fclose(fileID);

% This should be the right order (checked)
tt = reshape(gridtot,[gridparams{3} gridparams{2} gridparams{1}]);
tt = permute(tt,[3 2 1]);

hdrgrid = [double(gridparams{1}) double(gridparams{2}) double(gridparams{3}) ...
    double(gridparams{4}) double(gridparams{5}) double(gridparams{6}) ...
    double(gridparams{7}) double(gridparams{8}) double(gridparams{9})];

grid_ori = [gridori{4} gridori{6}];
