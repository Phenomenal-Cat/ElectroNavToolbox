function ENT_NormalizeToNative(TemplateFile, NativeFile, AtlasFiles, MaskFile)

%========================= ENT_NormalizeToNative.m ========================
% This function runs SPM8's spatial normalization method [1], to warp the
% specified atlas volumes to an individual animal's brain in native space.
%
% INPUTS:
%   TemplateFile:   full path of T1-weighted template volume being warped
%   NativeFile:     full path of individual's volume being warped to
%   AtlasFiles:     atlas volume or array of atlas structure volumes to warp
%   MaskFile:       (optional) volume used to weight normalization
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
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

if exist('spm','file')~=2
    error('SPM not detected on MATLAB path!');
end
% [ENTroot, temp, ext] = fileparts(mfilename('fullpath'));
% addpath(genpath(ENTroot));

ENTroot = '/Volumes/PROJECTS/murphya/EN_data';

%============== TEMPORARY!!! Hardcoded file paths
TemplateFile = fullfile(ENTroot, 'Atlases/inia19/inia19-T1-brain.nii');
NativeFile = fullfile(ENTroot,'Subjects/Layla/Layla_GridScan_ACPC_BET_Masked.nii');
% AtlasFiles = fullfile(ENTroot, 'Atlases/inia19/inia19-NeuroMaps.nii');
MaskFile = fullfile(ENTroot,'Subjects/Layla/Layla_ACPC_LesionMask.nii');

%============== Check inputs
if ~exist('TemplateFile','var') || ~exist(TemplateFile,'file')
    [file,path] = uigetfile('*.nii;*.img','Select template volume to warp from');
    TemplateFile = fullfile(path, file);
end
if ~exist('NativeFile','var') || ~exist(NativeFile,'file')
    [file,path] = uigetfile('*.nii;*.img','Select native volume to warp to');
    NativeFile = fullfile(path, file);
end
if ~exist('AtlasFiles','var') || ~exist(AtlasFiles,'file')
    AtlasFiles = uipickfiles('REFilter','.nii;.img', 'Prompt', 'Select volumes to transform');
end
if ischar(AtlasFiles)
    AtlasFiles = {AtlasFiles};
end
    
%=============== Calculate paramaters of native space
NativeNii = load_nii(NativeFile);                                               % load native volume
NativeNii.VoxDim = NativeNii.hdr.dime.pixdim(2:4);                              % get pixel dimensions (mm)
NativeNii.Origin = NativeNii.hdr.hist.originator(1:3);                          % get origin
NativeNii.BB(1,:) = -NativeNii.Origin.*NativeNii.VoxDim;                        % Calculate bounding box
NativeNii.BB(2,:) = (size(NativeNii.img)-NativeNii.Origin).*NativeNii.VoxDim;   

%=============== Estimate non-linear transform
NativeS = spm_vol(NativeFile);                                                  % Get structure for Native MRI volume
NativeH = spm_read_vols(NativeS);                                               % Load Native MRI volume
TemplateH = spm_read_vols(spm_vol(TemplateFile));
Matname = fullfile(ENTroot,'Subjects/Layla/INIA19_to_Layla_sn.mat');  
if exist('MaskFile','var')
    NativeMaskH = spm_read_vols(spm_vol(MaskFile));
else
    NativeMaskH = [];
end  
flags.smosrc = 0;                   % template image is already smoothed, so don't smooth again.
flags.smosrcsmoref = 4;             % smooth native volume to match template (4 mm FWHM Gaussian).
flags.smosrcregtype = 'rigid';      % regularisation type for affine registration. Use 'none' or 'rigid'
flags.smosrccutoff = 25;           	% Cutoff of the DCT bases.  Lower values mean more basis functions are used
flags.smosrcnits = 16;              % number of nonlinear iterations
flags.smosrcreg = 1;                % amount of regularisation. Default = 1. Increase if distortions appear
spm_normalise(NativeH, TemplateH, Matname, NativeMaskH, [], flags);

%=============== Apply estimated transform to images
flags.preserve = 0;
flags.bb = NativeNii.BB;            % Destination bounding box (mm)
flags.vox = NativeNii.VoxDim;     	% Destination voxel size (mm)
flags.Interp = 1;                   % Interpolation, 1 = trilinear?
flags.Wrap = [0 0 0];               % don't apply wrapping
flags.prefix = 'warped_';           % prefix added to warped volumes

for v = 1:numel(AtlasFiles)
    VO = spm_write_sn(AtlasFiles{v}, Matname, flags);
    V1 = spm_write_vol(VO, VO.dat);
end