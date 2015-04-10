
%========================= NormalizeToNative ==============================
% This function runs SPM8's spatial normalization method [1], to warp the
% specified atlas to an individual animal's brain in native space.
%
% DEFAULT PARAMATERS:
%   NativeT1 smoothing:     2-4mm FWHM (based on McLaren et al., 2010)
%   
%
% REFERNCES:
%   Ashburner, J., & Friston, K. J. (1999). Nonlinear spatial normalization 
%       using basis functions. Human brain mapping, 7(4), 254-266.
%   McLaren DG, Kosmatka KJ, Kastman EK, Bendlin BB, Johnson SC(2010). Rhesus 
%       macaque brain morphometry: A methodological comparison of voxel-wise 
%       approaches. Methods, 50(3): 157?165.
%   Yassa MA & Stark CEL (2009). A quantitative evaluation of cross-participant 
%       registration techniques for MRI studies of the medial temporal
%       lobe. NeuroImage, 44(2): 319?327.
%   http://www.fil.ion.ucl.ac.uk/spm/doc/books/hbf2/pdfs/Ch3.pdf
%   https://github.com/ritcheym/fmri_misc/blob/master/batch_template.m
%
% REVISIONS:
%   29/08/2014 - Written by APM
%==========================================================================

if exist('spm','file')~=2
    error('SPM not detected on MATLAB path!');
end
[ENTroot, temp, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(ENTroot));

%============== TEMPORARY!!! Hardcoded file paths 
AtlasFile = fullfile(ENTroot, 'Atlases/inia19/inia19-NeuroMaps.nii');
TemplateFile = fullfile(ENTroot, 'Atlases/inia19/inia19-T1-brain.nii');
NativeFile = fullfile(ENTroot,'Subjects/Layla/Layla_GridScan_ACPC_BET_Masked.nii');
MaskFile = fullfile(ENTroot,'Subjects/Layla/Layla_ACPC_LesionMask.nii');

%=============== Calculate paramaters of native space
NativeNii = load_nii(NativeFile);
NativeNii.VoxDim = NativeNii.hdr.dime.pixdim(2:4);    % get pixel dimensions
NativeNii.Origin = NativeNii.hdr.hist.originator(1:3);
NativeNii.BB(1,:) = -NativeNii.Origin.*NativeNii.VoxDim;
NativeNii.BB(2,:) = (size(NativeNii.img)-NativeNii.Origin).*NativeNii.VoxDim;

%=============== Estimate non-linear transform
NativeS = spm_vol(NativeFile);                      % Get structure for Native MRI volume
NativeH = spm_read_vols(NativeS);                   % Load Native MRI volume
NativeMaskH = spm_read_vols(spm_vol(MaskFile));           
TemplateH = spm_read_vols(spm_vol(TemplateFile));
Matname = fullfile(ENTroot,'Subjects/Layla/INIA19_to_Layla_sn.mat');

flags.smosrc = 0;                   % template image is already smoothed, so don't smooth again.
flags.smosrcsmoref = 8;             % smooth native volume to match template (8 mm FWHM Gaussian).
flags.smosrcregtype = 'rigid';      % regularisation type for affine registration. Use 'none' or 'rigid'
flags.smosrccutoff = 25;           	% Cutoff of the DCT bases.  Lower values mean more basis functions are used
flags.smosrcnits = 16;              % number of nonlinear iterations
flags.smosrcreg = 1;                % amount of regularisation. Default = 1. Increase if distortions appear
spm_normalize(NativeH, TemplateH, Matname, NativeMaskH, [], flags);

%=============== Apply estimated transform to images

flags.preserve = 0;
flags.bb = NativeNii.BB;            % Destination bounding box (mm)
flas.vox = NativeNii.VoxDim;        % Destination voxel size (mm)
flags.Interp = 1;                   % Interpolation, 1 = trilinear?
flags.Wrap = [0 0 0];               % don't apply warpping
flags.prefix = 'warped_';           
