% Filenames = CreateStructureMask(StructIndex, StructNames, AtlasFile, OutputDir)

%===================== ENT_CreateStructureMask.m ==========================
% Isolates the structure(s) specified by the indices 'StructIndex' from an
% atlas volume (default = NeuroMaps) and saves them as a binary mask.
%
% INPUTS:
%   StructIndex:    A cell containing structure indices (1 or 2 per structure)
%   StructNames:    A cell containing file names for each structure
%   AtlasFile:      full path of atlas volume to extract structures from.
%   OutputDir:    	full path to save new volumes to.
%   
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function Filenames = ENT_CreateStructureMask(StructIndex, StructNames, AtlasFile, OutputDir)

if nargin == 0
    [StructIndex, StructNames] = ENT_GetStructureIndex;
    if isempty(StructIndex)
        return;
    end
end
if exist('OutputDir','var')==0
    OutputDir = cd;
end
if exist('AtlasFile','var')==0 || exist(AtlasFile,'file')==0
    [file, path] = uigetfile('*.nii;*.img', 'Select atlas volume file');
    AtlasFile = fullfile(path, file);
end
AtlasNii = load_nii(AtlasFile);                                     % Load nifti file
for s = 1:numel(StructIndex)                                        % For each structure index provided...
    StructNii = AtlasNii;                                           % Copy nii struct
    StructNii.img = double(StructNii.img);                          % Convert volume class to doubles
    StructNii.img(~ismember(StructNii.img, StructIndex{s})) = 0;    % Set all voxels outside selected structure to zero
    StructNii.img(ismember(StructNii.img, StructIndex{s})) = 1;     % Set all voxels inside selected structure to one
    Filenames{s} = fullfile(OutputDir,[StructNames{s},'.nii']);     
    save_nii(StructNii, Filenames{s});                              % Save new structure mask as nifti volume
    fprintf('Structure volume %s.nii has been saved!\n', StructNames{s});
end

end