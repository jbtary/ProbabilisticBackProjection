function [ ] = plot_bkproj_prob_2d( rang, prod_prob, sta, title_name, hidplot )
% =========================================================================
% Function: plot_bkproj_prob_2d.m
%
% Plot results of back projected probabilities
% Auto choose slices passing through max of proba distribution
% =========================================================================
% Define size of station marker
sta_size = 7;

% Get max of probability cube
[~,I] = max(prod_prob(:));
[iy,ix] = ind2sub(size(prod_prob),I); clear I

maxp = [(ix-1)*rang(1,3) (iy-1)*rang(2,3)]; % Max proba position

if hidplot == 1; set(gcf,'Visible','Off'); else 
    figure; end

if max(prod_prob(:)) > 1e3
    imagesc((rang(2,1):rang(2,3):rang(2,2)),(rang(1,1):rang(1,3):rang(1,2)),log10(prod_prob))
else
    imagesc((rang(2,1):rang(2,3):rang(2,2)),(rang(1,1):rang(1,3):rang(1,2)),prod_prob)
end

hold on
plot(sta(:,1),sta(:,2),'v','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',sta_size)
plot(maxp(1),maxp(2),'w*','MarkerSize',10)
plot(maxp(1),maxp(2),'ko','MarkerSize',10)
ylabel('y (km)'); xlabel('x (km)')
title(title_name,'Interpreter','None')
colormap(jet)
set(gca,'YDir','normal')
colorbar


