% See events containing only hydrophone data
% Inputs:
% path: path to the directory where you want to start looking for the event
% fmin: minimum frequency in Hz for filtering (butterworth filter order 6)
% fmax: maximum frequency in Hz for filtering

function seeVAevent(path,fmin,fmax)

[file,path2] = uigetfile(path);
load([path2 file])

%%%%%
% Plot picks if there are some
if isfield(hdr,'ppicks')
    for jj=1:size(hdr.stap,1)
        compp{jj}=hdr.stap(jj,:); ppicks(jj)=hdr.ppicks(jj);
    end
    for jj=1:size(hdr.stas,1)
        comps{jj}=hdr.stas(jj,:); spicks(jj)=hdr.spicks(jj);
    end
end
%%%%%

figure(35); hold on; set(gcf,'Position', [10 10 1200 800])

maxtime = 0;
for jj=1:length(dath)
    if length(hdr.st) > 1
        fs = 1/hdr.st(jj);
    else
        fs = 1/hdr.st;
    end
    
    [b,a] = butter(6,[fmin fmax]/(fs/2),'bandpass');
    datfh = filtfilt(b,a,dath{jj,1});
    time = (0:length(datfh)-1)*(1/fs);
    plot(time,datfh/max(abs(datfh))+jj,'k');
    maxtmp = max(time);
    if maxtmp > maxtime; maxtime = maxtmp; end
    
    %%%%%
    % Plot the picks
    if strcmp(hdr.array{jj},'SAP') == 1
        comp{jj} = [hdr.array{jj}(1:2) hdr.comp{jj}];
    else
        comp{jj} = [hdr.array{jj} hdr.comp{jj}];
    end

    if exist('ppicks','var') == 1; [~,~,ia] = intersect(comp{jj},compp);
    else ia = []; end
    if exist('spicks','var') == 1; [~,~,ib] = intersect(comp{jj},comps);
    else ib = []; end

    if isempty(ia) ~= 1
        pick = (ppicks(ia) - hdr.tbeg(jj))*86400;
        plot(pick,jj,'ko','MarkerFaceColor','y'); clear pick
    end
    if isempty(ib) ~= 1
        pick = (spicks(ib) - hdr.tbeg(jj))*86400;
        plot(pick,jj,'ko','MarkerFaceColor','g'); clear pick
    end
    clear ia ib pick comp
    %%%%%

    clear datf* time maxtmp
end

axis([0 maxtime 0 jj+1]); box on;
xlabel('Time (s)'); ylabel('Array-sensor numbers')
title(['File ' file],'Interpreter','none');
set(gca,'YTick',1:jj);

% Set axis with instrument's names
for jj=1:length(hdr.comp); comp{jj} = [hdr.array{1,jj} '-' hdr.comp{1,jj}]; end
set(gca,'YTickLabel',comp)
