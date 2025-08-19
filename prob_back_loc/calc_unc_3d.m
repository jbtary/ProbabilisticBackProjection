% Calculation of uncertainties at a given level of the max of the PDF
% 3D
% 
% Input:
%   prod_prob: location PDF
%   fact: factor by which to divide the max of exp(prod_prob) to define the
%   uncertainty zone (needs to be > 1, ex: 2)
%   rang: grid variable with (x0 xmax dx, y0 ymax dy, z0 zmax dz)
% 
% Output:
%   uncs: uncertainties in x, y, and z (same units as rang)
%   uncs_coord: coordinates of the uncertainty zone (same units as rang)

function [uncs,uncs_coord] = calc_unc_3d(prod_prob,fact,rang)

prod_log = prod_prob - (min(prod_prob(:)));
prod_log(isnan(prod_log)) = 0;

% Transform the PDF into a vector, keep track of indexes
pdf_vec = prod_log(:);

% Sort the vector in descending order, keep track of indexes
[pdf_sort,inds] = sort(pdf_vec,'descend');

% % Calculate the cum sum of the vector
% pdf_cum = cumsum(pdf_sort);
% % Cut the cum sum at the given percentage
% iv = find(pdf_cum > perc*pdf_cum(end),1);

% Get all values above a given value
iv = find(pdf_sort > max(pdf_sort)-log(fact),1,'last');

% Return the uncertainties
pdf_fin = zeros(size(prod_prob));
pdf_fin = pdf_fin(:);
pdf_fin(inds(1:iv)) = 1;
pdf_fin = reshape(pdf_fin,[size(prod_prob,1) size(prod_prob,2) size(prod_prob,3)]);

[iz,iy,ix] = ind2sub(size(pdf_fin),find(pdf_fin > 0));
uncs = [(max(ix)-min(ix))*rang(1,3) (max(iy)-min(iy))*rang(2,3) (max(iz)-min(iz))*rang(3,3)];
uncs_coord = [(min(ix)-1)*rang(1,3) (max(ix)-1)*rang(1,3);...
    (min(iy)-1)*rang(2,3) (max(iy)-1)*rang(2,3); (min(iz)-1)*rang(3,3) (max(iz)-1)*rang(3,3)];
