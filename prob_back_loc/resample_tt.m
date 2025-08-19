% Function to resample the matrix of travel-times
% Better to have grids that can easily be combined (compatible spacings,
% same origin, new grid needs to be contained in previous grid)

function sam_resam = resample_tt(sam,hdrgrid,dx,maxx,maxy,maxz)

% We will assume that either all intervals are larger than previous grid or
% that they are all greater. Taking into account all cases would be too much
% for my case.
if dx > hdrgrid(7) % Preferred case
    % If the grid spacings are incompatible and need interpolation 
    % (not taken into account either)
    if mod(dx,hdrgrid(7)) ~= 0
        disp('Incompatible grid spacings')
        sam_resam = [];
    else
        rat = dx/hdrgrid(7);
        % Resample as if we were using nearest neighbour algo
        disp('Travel-time grid resampling')
        sam_resam = sam(1:rat:(maxx/hdrgrid(7))+1,...
            1:rat:(maxy/hdrgrid(8))+1,1:rat:(maxz/hdrgrid(9))+1);
    end
else
    disp('Come on, use a better travel-time grid! Resampling anyway.')
    % Needs Image Proc toolbox
    rat = hdrgrid(7)/dx;
    sam_resam = imresize3(sam,rat,'linear');
end
