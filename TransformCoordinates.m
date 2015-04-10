function xyz2 = TransformCoordinates(xyz, WarpFile, Inv)

%====================== TransformCoordinates.m ============================
%
%
% INPUTS:
%   xyz:        n x 3 matrix of coordinates
%   WarpFile:   full filename of .mat file produced by SPM8 normalization
%               of volume A to volume B (typically with the suffix '_sn.mat').
%   Inv:        Optional flag. Default = 0. 1 = perform inverse transformation.
%
% EXAMPLE:
%   xyz2 = TransformCoordinates(xyz, WarpFile, 0);
%
% REVISIONS:
%   12/10/2014 - Written by Aidan Murphy (murphyap@mail.nih.gov) based on 
%                get_orig_coord2.m by John Ashburner.
%==========================================================================

spm_check_installation;
if size(xyz,2)~=3, error('xyz must be an N x 3 matrix'); end;
xyz = xyz';

load(WarpFile);
if exist('Inv','var')
    if Inv == 1
        Mat = inv(VG.mat);
    else
        Mat = VG.mat;
    end
else
    Mat = VG.mat;
end
xyz = Mat(1:3,:)*[xyz ; ones(1,size(xyz,2))];
d   = VG.dim(1:3);
Mult = VF.mat*Affine;

if (prod(size(Tr)) == 0),
    affine_only = 1;
    basX = 0; tx = 0;
    basY = 0; ty = 0;
    basZ = 0; tz = 0;
else
    affine_only = 0;
    basX = spm_dctmtx(d(1),size(Tr,1),xyz(1,:)-1);
    basY = spm_dctmtx(d(2),size(Tr,2),xyz(2,:)-1);
    basZ = spm_dctmtx(d(3),size(Tr,3),xyz(3,:)-1);
end

if affine_only,
	xyz2 = Mult(1:3,:)*[xyz ; ones(1,size(xyz,2))];
else
	for i=1:size(xyz,2)
		bx = basX(i,:);
		by = basY(i,:);
		bz = basZ(i,:);
		tx = reshape(...
			reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3))...
			*bz', size(Tr,1), size(Tr,2) );
		ty = reshape(...
			reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3))...
			*bz', size(Tr,1), size(Tr,2) );
		tz =  reshape(...
			reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3))...
			*bz', size(Tr,1), size(Tr,2) );
		xyz2(:,i) = Mult(1:3,:)*[xyz(:,i) + [bx*tx*by' ; bx*ty*by' ; bx*tz*by']; 1];
    end
end
xyz2 = xyz2';
