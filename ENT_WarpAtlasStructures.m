function ENT_WarpAtlasStructures(AtlasFile, Transform)

%====================== ENT_WarpAtlasStructures.m =========================
% This wrapper function allows the user to select specific anatomical
% structures from a specified atlas and create volumetric and surface mesh
% representations of them that are spatially normalized to an individual
% monkey ('native space'). 
%
% INPUTS:
%   AtlasFile: 	full path of the atlas volume (.nii) file
%   Transform:  full path of the transform matrix (.mat) created by 
%               spm_normalize.m for non-linearly warping the atlas to an
%               individual animal's native (e.g. AC-PC) space.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

%============ CHECK INPUTS
if ~exist('AtlasFile', 'var')
    [file, path] = uigetfile('*.nii;*.img', 'Select atlas volume file');
    AtlasFile = fullfile(path, file);
end
if ~exist('Transform', 'var')
    [file, path] = uigetfile('*.mat', 'Select transform matrix file');
    Transform = fullfile(path, file);
end

%============ SELECT BRAIN STRUCTURES
Bilateral = 1;      
[StructIndex, StructNames] = ENT_GetStructureIndex([], Bilateral);          % Get atlas indices for requested structure(s)

if numel(StructNames)>1
    StructIndex{end+1} = cell2mat(StructIndex);                            	% Combine all structures?
    StructNames{end+1} = ['All_',Structures{1}];                          	% Give combined name
end

%============ CREATE OUTPUT DIRECTORY
SubjectDir = uigetdir('','Select subject directory');
StructVolDir = fullfile(SubjectDir, 'StructureVolumes');
if exist(StructVolDir,'dir')==0
    mkdir(StructVolDir);
end                                                     
StructFilenames = ENT_CreateStructureMask(StructIndex, StructNames, AtlasFile, StructVolDir); 	% Create the structure volumes
StructFilenames = ENT_ApplyTransform(StructFilenames, Transform);           % Transform atlas structures to native (AC-PC) space
ENT_VolumeToSurface(StructFilenames, 'vtk', OutputDir, 1);                  % Create surface meshes of structures, view in 3D and save