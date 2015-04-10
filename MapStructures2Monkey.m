
%========================= MapStructureToMonkey.m =========================
% This wrapper function allows the user to select specific anatomical
% structures from a specified atlas and create volumetric and surface mesh
% representations of them that are spatially normalized to an individual
% monkey ('native space'). 
%
% REQUIREMENTS:
%   ElectroNavToolbox:  https://github.com/MonkeyGone2Heaven/ElectroNavToolbox
%   SPM8:               http://www.fil.ion.ucl.ac.uk/spm/software/spm8/
%==========================================================================


Structures = {'ventral_lateral'};
SubjectDir = '/Volumes/PROJECTS-1/murphya/Toolboxes/ElectroNavToolbox/Subjects/Layla';
SubjectTransform = 'inia19-NeuroMaps_sn.mat';


Bilateral = 1;      
[StructIndex, StructNames] = GetStructureIndex(Structures, Bilateral);      % Get atlas indices for requested structure(s)

if numel(StructNames)>1
    StructIndex{end+1} = cell2mat(StructIndex);                            	% Combine all structures?
    StructNames{end+1} = ['All_',Structures{1}];                          	% Give combined name
end

%============ Create directory to save nifti volumes of anatomical structures
StructVolDir = fullfile(SubjectDir, 'StructureVolumes');
if exist(StructVolDir,'dir')==0
    mkdir(StructVolDir);
end                                                     
Filenames = CreateStructureMask(StructIndex, StructNames, StructVolDir);              % Create the structure volumes


%=========== Transform atlas structures to native ACPC



%=========== Create .VTKs of structures and view in 3D
fh = RenderVolumes(StructVolDir);
