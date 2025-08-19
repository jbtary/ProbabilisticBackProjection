function [ ] = plot_bkproj_prob_3d( rang, prod_prob, sta, title_name, hidplot )
% =========================================================================
% Function: plot_bkproj_prob_3d.m
%
% Plot results of back projected probabilities
% Auto choose slices passing through max of proba distribution
% =========================================================================
% Define size of station marker
sta_size = 7;

% Get max of probability cube
[~,I] = max(prod_prob(:));
[iz,iy,ix] = ind2sub(size(prod_prob),I); clear I

% Get 2D arrays
sliceXZ = squeeze(prod_prob(:,iy,:)); % Cross-section with constant Y
sliceYZ = squeeze(prod_prob(:,:,ix)); % Cross-section with constant X
sliceXY = squeeze(prod_prob(iz,:,:)); % Map
maxp = [(ix-1)*rang(1,3) (iy-1)*rang(2,3) (iz-1)*rang(3,3)]; % Max proba position

if hidplot == 1; set(gcf,'Visible','Off'); else 
    figure; end

set(gcf,'Position', [10 10 800 600]) % Resize figure size to make it bigger
axes('position',[.1 .5 .5 .4]); % Pos minx miny sizex sizey, Hz slice
if max(prod_prob(:)) > 1e3
    imagesc(rang(1,1):rang(1,3):rang(1,2),rang(2,1):rang(2,3):rang(2,2),log10(sliceXY))
else
    imagesc(rang(1,1):rang(1,3):rang(1,2),rang(2,1):rang(2,3):rang(2,2),sliceXY)
    clim([max(prod_prob(:))-log(10) max(prod_prob(:))])
end
hold on
plot(sta(:,1),sta(:,2),'v','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',sta_size)
plot(maxp(1),maxp(2),'w*','MarkerSize',10)
plot(maxp(1),maxp(2),'ko','MarkerSize',10)
ylabel('y (km)')
title(title_name,'Interpreter','None','BackgroundColor', [.9 .9 .9])
c_map = colormap(jet); % Set background color as white
c_map(1,:) = [1 1 1];
colormap(c_map)
set(gca,'YDir','normal'); set(gca,'FontSize',12)

axes('position',[.1 .3 .5 .15]); % West-East cross-section
if max(prod_prob(:)) > 1e3
    imagesc((rang(1,1):rang(1,3):rang(1,2)),(rang(3,1):rang(3,3):rang(3,2)),log10(sliceXZ))
else
    imagesc((rang(1,1):rang(1,3):rang(1,2)),(rang(3,1):rang(3,3):rang(3,2)),sliceXZ)
    clim([max(prod_prob(:))-log(10) max(prod_prob(:))])
end
hold on
plot(maxp(1),maxp(3),'w*','MarkerSize',10)
plot(maxp(1),maxp(3),'ko','MarkerSize',10)
xlabel('x (km)'); ylabel('z (km)'); set(gca,'FontSize',12)

axes('position',[.65 .5 .15 .4]); % North-South cross-section
if max(prod_prob(:)) > 1e3
    imagesc((rang(3,1):rang(3,3):rang(3,2)),(rang(2,1):rang(2,3):rang(2,2)),log10(sliceYZ'))
else
    imagesc((rang(3,1):rang(3,3):rang(3,2)),(rang(2,1):rang(2,3):rang(2,2)),sliceYZ')
    clim([max(prod_prob(:))-log(10) max(prod_prob(:))])
end
hold on
plot(maxp(3),maxp(2),'w*','MarkerSize',10)
plot(maxp(3),maxp(2),'ko','MarkerSize',10)
xlabel('z (km)')
set(gca,'YDir','normal'); set(gca,'FontSize',12)
colorbar
