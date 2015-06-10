
%========================== ENT_ApplyTransform.m ==============================
% This function applies the spatial transformation specified by a 4x4 
% matrix and applies it to the specified volume(s).
%
% INPUTS:
%   Volumes:        full path of the MR volume(s) to apply transform to.
%   TargetVolume:   full path of target volume (native animal space)
%   Transform:      full path of .mat file containing the 4x4 transform matrix 
%                   created by spm_normalize.m
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function ENT_ApplyTransform(Volumes, TargetVolume, Transform)

%========== CHECK REQUIRED FUNCTIONS
if ~exist('spm.m','file')
    error(  ['ERROR: SPM was not detected on the MATLAB path! ',...
             'ENT_ApplyTransform.m requires SPM''s ''spm_write_sn.m'' function.']);
end

%========== CHECK INPUTS
if ~exist('Volumes','var')
    Volumes = uipickfiles('REFilter','.nii;.img', 'Prompt', 'Select volumes to transform');
end
if ischar(Volumes)
    Volumes = {Volumes};
end
if ~exist('TargetVolumes','var')
    [file, path] = uigetfile('*.nii;*.img', 'Select target volume');
    TargetVolume = fullfile(path, file);
end
if ~exist('Transform','var') || ~exist(Transform,'file')
    [file, path] = uigetfile('*.mat', 'Select transform matrix file');
    Transform = fullfile(path, file);
end
% if ischar(Transform)
%     temp = load(Transform);
%     T = temp.Affine;
% end

%========== APPLY TRANSFORM AND SAVE
VI = spm_vol(TargetVolume);                   	% Load target volume
Flags.interp   	= 1;                            % Interp method = 
Flags.wrap      = [];                           % No wrapping
Flags.vox       = diag(VI.mat(1:3,1:3))';       % Get voxel dimensions
Flags.bb        = [(-VI.dim.*Flags.vox)/2;    	% Set bounding box (mm)
                   (VI.dim.*Flags.vox)/2];
Flags.preserve	= 0;                            % No volume correction
    
for V = 1:numel(Volumes)
    VO = spm_write_sn(Volumes{V}, Transform, Flags);
    V1 = spm_write_vol(VO, VO.dat);
end