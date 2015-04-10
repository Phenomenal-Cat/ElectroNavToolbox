
%==========================================================================
% Apply a binary mask volume to another volume of identical size. 
%
% INPUTS:
%   MaskNiiFile:    Full path of binary mask .nii file 
%   VolNiiFile:     Full path of volume .nii file
%   Subtract:       0 = all voxels of intensity zero in the mask image
%                   become zero in the volume image.
%                   1 = all voxels of intensity one in the mask image
%                   become zero in teh volume image.
%
%==========================================================================

function ApplyMask(MaskNiiFile, VolNiiFile, Subtract)

if nargin == 0
    [Filename,Pathname,Indx] = uigetfile('*.nii','Select mask file');
    MaskNiiFile = fullfile(Pathname, Filename);
    [Filename,Pathname,Indx] = uigetfile('*.nii','Select volume to mask');
    VolNiiFile = fullfile(Pathname, Filename);
    Subtract = 1;
end
MaskNii = load_nii(MaskNiiFile);
VolNii = load_nii(VolNiiFile);
if size(MaskNii.img)~= size(VolNii.img)
    error('Mask and volume to mask are not the same size!\n');
end
if Subtract == 1
    VolNii.img(MaskNii.img == 1) = 0; 
elseif Subtract == 0
    VolNii.img(MaskNii.img == 0) = 0; 
end
VolNiiFilename = [VolNiiFile(1:end-4),'_Masked.nii'];
save_nii(VolNii, VolNiiFilename);