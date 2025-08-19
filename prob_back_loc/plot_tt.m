% Simple function to plot the travel-time tables
% nx,ny,nz: number of grid nodes
% dx,dy,dz: spacing between grid nodes
% ix,iy,iz: index of the nodes for sections (lower than number of grid nodes)
% sta: station positions in local coordinates (x,y,z)

function plot_tt(tt,nx,ny,nz,dx,dy,dz,ix,iy,iz,sta)

% Get 2D arrays
sliceXZ = squeeze(tt(:,iy,:)); % Cross-section with constant Y
sliceYZ = squeeze(tt(ix,:,:)); % Cross-section with constant X
sliceXY = squeeze(tt(:,:,iz)); % Map

figure
axes('position',[.1 .5 .5 .4]); % Pos minx miny sizex sizey, Hz slice
imagesc(0:dx:(nx-1)*dx,0:dy:(ny-1)*dy,sliceXY')
hold on
plot(sta(:,1),sta(:,2),'v','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',7)
ylabel('y (km)')
title('Travel-times in sec.')
colormap(jet)
set(gca,'YDir','normal')

axes('position',[.1 .15 .5 .3]); % West-East cross-section
imagesc(0:dx:(ny-1)*dx,0:dz:(nz-1)*dz,sliceXZ')
xlabel('x (km)'); ylabel('z (km)')

axes('position',[.65 .5 .25 .4]); % North-South cross-section
imagesc(0:dz:(nz-1)*dz,0:dy:(ny-1)*dy,sliceYZ)
xlabel('z (km)')
set(gca,'YDir','normal')
colorbar
