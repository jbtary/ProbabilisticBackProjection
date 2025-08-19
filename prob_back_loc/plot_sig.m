% Plotting function for signals
% 
% Input:
%   sig: signals to plot (Nsta x N samples)
%   t_axis_seg: time in sec. (Nsta x N samples)
%   sta_name: names of stations
%   ind: indexes of the stations to plot

function plot_sig(sig,t_axis_seg,sta_name,ind)

figure; hold on
ct = 1;
for ii = ind
    plot(t_axis_seg(ii,:),sig(ii,:)/max(sig(ii,:)) + (ct-1),'k')
    text(t_axis_seg(ii,1),ct-0.5,num2str(max(sig(ii,:))))
    ct = ct + 1;
end

ct = 1;
for ii = ind
    staplot{ct} = sta_name{ii}; ct = ct + 1;
end
set(gca,'YTick',0:ct-2)
set(gca,'YTickLabel',staplot)
ylim([-1 length(ind)])
xlabel('Time (s)')
ylabel('Station names')
title('Signals')
box on; grid on; set(gca,'FontSize',12)
