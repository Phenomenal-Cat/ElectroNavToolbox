% Filenames = ENT_CreateStructureMask(StructIndex, StructNames, AtlasName, OutputDir)

%===================== ENT_CreateStructureMask.m ==========================
% Isolates the structure(s) specified by the indices 'StructIndex' from an
% atlas volume and saves them as a new volume containing a binary mask.
%
% INPUTS:
%   StructIndex:    A cell containing structure indices (1 or 2 per structure)
%   StructNames:    A cell containing file names for each structure
%   AtlasName:      name of atlas to extract structures from.
%   OutputDir:    	full path to save new volumes to.
%   
% EXAMPLE:
%   OutputDir = '/Volumes/PROJECTS/murphya/EN_data/Atlases/D99/VTKs';
%   Filenames = ENT_CreateStructureMask({[196, 197, 198]}, 'D99_pulvinar_ROI', 'Saleem', OutputDir)
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function Filenames = ENT_CreateStructureMask(StructIndex, StructNames, AtlasName, OutputDir)

if nargin == 0
    Struct      = ENT_GetStructureIndex;
    StructIndex = {Struct.Indices};
    StructNames = {Struct.Names};
    if isempty(StructIndex)
        return;
    end
end
if exist('OutputDir','var')==0
    OutputDir = uigetdir(cd, 'Select path to save new volumes to');
end

addpath(genpath('/Volumes/PROJECTS/murphya/EN_data/EN_Atlases'));

%============= LOAD ATLAS DATA
Atlas = ENT_LoadAtlasData;
if ~exist('AtlasName','var')
    Selection   = listdlg('ListString',{Atlas.name},'SelectionMode','single','Listsize',[150, 100],'PromptString','Select atlas');
    AtlasName   = Atlas(Selection).name; 
else
    Selection   = find(~cellfun(@isempty, strfind(lower({Atlas.name}), lower(AtlasName)))); 
end
AtlasFile       = Atlas(Selection).AtlasFile;
AtlasStruct     = ENT_GetStructureIndex(AtlasName);

if ~exist(AtlasFile,'file')
    [file, path]= uigetfile('*.nii;*.img', 'Select atlas volume file');
    AtlasFile   = fullfile(path, file);
end
AtlasNii        = load_nii(AtlasFile);                                      % Load nifti file
if any(~isinteger(AtlasNii.img(:)))
    AtlasNii.img = round(AtlasNii.img);
end

%============= CREATE NEW VOLUME FOR EACH STRUCTURE
if ~iscell(StructNames)
    StructNames = {StructNames};
end
for s = 1:numel(StructIndex)                                            	% For each structure index provided...
    StructNii       = AtlasNii;                                             % Copy nii struct
    StructNii.img   = double(StructNii.img);                                % Convert volume class to doubles
    StructNii.img(~ismember(AtlasNii.img, StructIndex{s}))  = 0;            % Set all voxels outside selected structure to zero
    StructNii.img(ismember(AtlasNii.img, StructIndex{s}))   = 1;            % Set all voxels inside selected structure to one
    Filenames{s}    = fullfile(OutputDir,[StructNames{s},'.nii']);     
    save_nii(StructNii, Filenames{s});                                      % Save new structure mask as nifti volume
    fprintf('Structure volume %s.nii has been saved!\n', StructNames{s});
end

end