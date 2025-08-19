% Function to calculate the max time difference between two stations
% (=direct ray path between them)
% 
% Here we assume that the 0,0 of the travel-time grid is the same as the
% 0,0 used to calculate local station coordinates (else one would need to
% adjust for this)
% 
% Input:
%   ttsta1: travel-time grid of station 1 in sec. (nx,ny,nz)
%   ttsta2: travel-time grid of station 2 in sec. (nx,ny,nz)
%   hdrgrid: grid paramaters (nx,ny,nz,xshift,yshift,zshift,dx,dy,dz)
%   vel: 3d velocity model with same grid parameters as travel-time grid
%   maxx,maxy,maxz: maximum grid dimensions for the location grid (assumed
%   same 0,0,0)
% 
% Output:
%   max_tlag: direct wave propagation time between 2 stations as calculated
%   by the eikonal equation in sec.
% 

function max_tlag = diff_time_tt(ttsta1,ttsta2,hdrgrid,vel,maxx,maxy,maxz)

nx = hdrgrid(1,1); ny = hdrgrid(1,2); nz = hdrgrid(1,3);
dx = hdrgrid(1,7); dy = hdrgrid(1,8); dz = hdrgrid(1,9);
maxgridx = (nx-1)*dx;
maxgridy = (ny-1)*dy;
maxgridz = (nz-1)*dz;

[indx,indy,indz] = size(ttsta1);
% Location grid needs to be equal or inside travel-time grids
if maxx < maxgridx || maxy < maxgridy || maxz < maxgridz
    % Use round to avoid incompatible grids, but this must be avoided and
    % will throw an error later
    indx = round(maxx/dx)+1;
    indy = round(maxy/dy)+1;
    indz = round(maxz/dz)+1;
end

if maxx > maxgridx || maxy > maxgridy || maxz > maxgridz
    disp('diff_time_tt.m, Location grid larger than travel-time grid, abort.')
    max_tlag = [];
    return
end

% Get max of the differential time on complete grid, add a small amount 
% to take into account the space between node and real position of station
diff_tt = abs(ttsta1(1:indx,1:indy,1:indz) - ttsta2(1:indx,1:indy,1:indz));
max_tlag = max(diff_tt(:)) + (sqrt(dx^2 + dy^2 + dz^2)/min(vel(vel>0.1)));
