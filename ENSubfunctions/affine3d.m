function [new_im] = affine3d(old_im, M, range_x, range_y, range_z, method)
%function [new_im] = affine3d(old_im, M, range_x, range_y, range_z, method)
% old_im: old image
% M: 4x4 homogeneous backward affine transformation matrix
% range_x: x-points in destination
% range_y: y-points in destination
% range_z: z-points in destination
% interp_method (optional): interpolation method'nearest', 'linear', 
% 'spline', 'cubic'. Add * for more speed, see help interp3 for details

% Author: Martijn Steenwijk
% Date: October 17, 2009

% % Usage example - simple translation
% load mri.mat
% D = squeeze(D(:,:,1,:));
% % Create affine transformation matrix, simply shift (x,y) = (50,25)
% M = [1 0 0 50; 0 1 0 25; 0 0 1 0; 0 0 0 1];
% % Invert M, since the interpolation is backward. Meanwhile subsample the 
% output volume x-direction by a factor of two.
% D_new = affine3d(D,inv(M),1:2:128,1:128,1:27);
% figure
% subplot(1,2,1)
% imagesc(D(:,:,10))
% title('Original volume')
% subplot(1,2,2)
% imagesc(D_new(:,:,10))
% title('Shifted volume')


if nargin < 5 || nargin > 6 
   error('Wrong number of input arguments');
elseif nargin == 5
    method = 'linear';
end

% convert old_im to single
old_im = single(old_im);

%Get all points in destination to sample
[yg xg zg] = meshgrid(range_y,range_x,range_z);
xyz = [reshape(xg,numel(xg),1)'; reshape(yg,numel(yg),1)'; reshape(zg,numel(zg),1)'];
xyz = [xyz; ones(1,size(xyz,2))];

%transform into source coordinates
uvw = M * xyz;

%Remove homogeneous
uvw = uvw(1:3,:)';

%Sample
xi = reshape(uvw(:,1), length(range_x),length(range_y),length(range_z));
yi = reshape(uvw(:,2), length(range_x),length(range_y),length(range_z));
zi = reshape(uvw(:,3), length(range_x),length(range_y),length(range_z));

% interp3 treats x and y in right-handed coordinate system, not in matrix
% index order, so we need to swap them here.
new_im = interp3(old_im,yi,xi,zi,method);

%Check for NaN background pixels - replace them with a background of 0
idx = find(isnan(new_im));
if(~isempty(idx))
    new_im(idx) = 0;
end
