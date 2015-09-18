function [dcm] = ENT_Rot2Mat(xrot, yrot, zrot)

% Takes roatation angles (degrees) about x, y and z axes and outputs a 3 x 3
% rotation matrix.

angles  = [xrot(:) yrot(:) zrot(:)];
dcm     = zeros(3,3,size(angles,1));
cang    = cos( angles );
sang    = sin( angles );

% Reverse (Z 1st, Y 2nd, X 3rd)
dcm(1,1,:) = cang(:,2).*cang(:,1);
dcm(1,2,:) = cang(:,2).*sang(:,1);
dcm(1,3,:) = -sang(:,2);
dcm(2,1,:) = sang(:,3).*sang(:,2).*cang(:,1) - cang(:,3).*sang(:,1);
dcm(2,2,:) = sang(:,3).*sang(:,2).*sang(:,1) + cang(:,3).*cang(:,1);
dcm(2,3,:) = sang(:,3).*cang(:,2);
dcm(3,1,:) = cang(:,3).*sang(:,2).*cang(:,1) + sang(:,3).*sang(:,1);
dcm(3,2,:) = cang(:,3).*sang(:,2).*sang(:,1) - sang(:,3).*cang(:,1);
dcm(3,3,:) = cang(:,3).*cang(:,2);


% % Reverse (X 1st, Y 2nd, Z 3rd)
% dcm(1,1,:) = cang(:,2).*cang(:,3);
% dcm(1,2,:) = sang(:,1).*sang(:,2).*cang(:,3) + cang(:,1).*sang(:,3);
% dcm(1,3,:) = -cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,3);
% dcm(2,1,:) = -cang(:,2).*sang(:,3);
% dcm(2,2,:) = -sang(:,1).*sang(:,2).*sang(:,3) + cang(:,1).*cang(:,3);
% dcm(2,3,:) = cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,3);        
% dcm(3,1,:) = sang(:,2);
% dcm(3,2,:) = -sang(:,1).*cang(:,2);
% dcm(3,3,:) = cang(:,1).*cang(:,2);

