% Plotting function for envelopes of xcorr
% 
% Input:
%   env: envelope of xcorr to plot
%   tlag1: lags corresponding to env
%   pair2: indexes of pairs of data corresponding to env
%   env_ids: indexes in env to select which ones to plot

function plot_env(env,tlag1,pair2,env_ids)

figure; hold on
ct = 1;
for ii = env_ids
    plot(tlag1,env(ii,:)/max(env(ii,:)) + (ct-1),'k')
    text(tlag1(1),ct-0.5,num2str(max(env(ii,:))))
    ct = ct + 1;
end

ct = 1;
for ii = env_ids
    pairs{ct} = num2str(pair2(ii,1:2)); ct = ct + 1;
end
set(gca,'YTick',0:ct-2)
set(gca,'YTickLabel',pairs)
xlabel('Lag (samples)')
ylabel('Indexes of pairs of stations')
title('Envelopes and max of envelopes')
box on; grid on; set(gca,'FontSize',12)
