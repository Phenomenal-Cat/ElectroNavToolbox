% Filenames = CreateStructureMask(StructIndex, StructNames, AtlasFile, Outdir)

%============================ GetStructureIndex ===========================
% Isolates the structure(s) specified by the indices 'StructIndex' from an
% atlas volume (default = NeuroMaps) and saves them as a binary mask.
%
% INPUTS:
%   StructIndex:    A cell containing structure indices (1 or 2 per structure)
%   StructNames:    A cell containing file names for each structure
%   AtlasFile:      full path of atlas volume to extract structures from.
%   Outdir:         full path to save new volumes to.
%   
% REVISION HISTORY:
%   2013 - Written by APM (murphyap@mail.nih.gov)
%==========================================================================

function Filenames = CreateStructureMask(StructIndex, StructNames, Outdir)

if nargin == 0
    [StructIndex, StructNames] = GetStructureIndex;
    if isempty(StructIndex)
        return;
    end
end
if exist('Outdir','var')==0
    Outdir = cd;
end
if exist('AtlasFile','var')==0
    AtlasFile = 'Atlases/inia19/inia19-NeuroMaps.nii';
end
if exist(AtlasFile,'file')==0
    error('Atlas file ''%s'' was not found!\n', AtlasFile);
end
AtlasNii = load_nii(AtlasFile);                                     % Load nifti file
for s = 1:numel(StructIndex)                                        % For each structure index provided...
    StructNii = AtlasNii;                                           % Copy nii struct
    StructNii.img = double(StructNii.img);                          % Convert volume class to doubles
    StructNii.img(~ismember(StructNii.img, StructIndex{s})) = 0;    % Set all voxels outside selected structure to zero
    StructNii.img(ismember(StructNii.img, StructIndex{s})) = 1;     % Set all voxels inside selected structure to one
    Filenames{s} = fullfile(Outdir,[StructNames{s},'.nii']);
    save_nii(StructNii, Filenames{s});                              % Save new structure mask as nifti volume
    fprintf('Structure volume %s.nii has been saved!\n', StructNames{s});
end

end