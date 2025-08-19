% Function to return station names, locations, and corresponding data, as
% well as order in master station file, based on which ones are present in 
% an event file or a sub-selection
% 
% Needs: deg2utm.m
% Input:
%   pathandfilename: path and filename of the event to load/locate (needs
%   to contain header with station info)
%   allsta: path and filename of the text file including all station info
%   (name and lat/long locations) of the survey that were used to calculate 
%   the mapping operator
%   grid0: lat long of the zero of the velocity grid (lower left corner)
%   (ex: [3 -84])
%   sta_sel: cell array of station names to select (can be empty if not
%   needed)
%   
% Output:
%   sta_name: names of the stations in the file to use
%   sta: location of the stations to use
%   ind: indexes of the stations to use in the complete list of stations
%   (to merge with mapping operator)
%   data: hydro data in same order as the complete list of stations
% 

function [sta_name,sta,ind,data,hdr2] = sta_extract_h(pathandfilename,allsta,grid0,sta_sel)

load(pathandfilename,'hdr','dath')

% Get station info: names and geo locations
for ii = 1:length(hdr.comp); comp{ii} = [hdr.array{ii} hdr.comp{ii}]; end
obslocs = importdata(allsta,' ');
stations = obslocs.textdata;
[sta_name,ind,icomp] = intersect(stations,comp,'stable');
sta_loc_geo = [obslocs.data(ind,2) obslocs.data(ind,1) obslocs.data(ind,3)]; % Lon Lat depth

% Convert geo to local coordinate system (utm) in km
[sta(:,1),sta(:,2)] = deg2utm(sta_loc_geo(:,2),sta_loc_geo(:,1));
sta = sta/1000; % Convert from m to km
[minx,miny] = deg2utm(grid0(1),grid0(2)); % Define the 0,0 (same as TT/Vel grid)
minx = minx/1000; miny = miny/1000; % Convert from m to km
sta(:,1) = sta(:,1) - minx;
sta(:,2) = sta(:,2) - miny;

if sum(sta_loc_geo(:,3)<0) > 0; sta(:,3) = -sta_loc_geo(:,3); end

% Reorder data in order of file with all station info
for ii = 1:length(sta_name)
    data.h{ii,1} = dath{icomp(ii),1};
    hdr2.tbeg(ii) = hdr.tbeg(icomp(ii));
end
% hdr2.hdr = hdr.hdr; 
hdr2.st = hdr.st(1);
if length(hdr.st) > 1
    for ii = 1:length(sta_name)
        hdr2.st(ii) = hdr.st(icomp(ii));
    end
end

% Only return a selection of stations if requested
if ~isempty(sta_sel)
    [sta_name,ia,ib] = intersect(sta_name,sta_sel,'stable');
    if length(sta_name) ~= length(sta_sel)
        disp('sta_extract.m, some stations selected were not present in event files')
    end
    sta = sta(ia,:);
    ind = ind(ia,:);
    datasv = data; clear data
    tbeg = hdr2.tbeg;
    hdr2 = rmfield(hdr2,'tbeg');
    for ii = 1:length(sta_name)
        data.h{ii,1} = datasv.h{ia(ii),1};
        hdr2.tbeg(ii) = tbeg(icomp(ii));
    end

    if length(hdr2.st) > 1
        st = hdr2.st;
        hdr2 = rmfield(hdr2,'st');
        for ii = 1:length(sta_name)
            hdr2.st(ii) = st(icomp(ii));
        end
    end
end
