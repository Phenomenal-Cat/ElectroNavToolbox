function ENT_MakeGridVolume(GridFile, HeaderSource, VoxelSize)

%========================== ENT_MakeGridVolume.m ==========================
% Converts a 3D model of a recording grid (saved in .stl/.obj/.ply/.vtk 
% formats), into a niftii (.nii) volume file. This is next manually aligned 
% to a structural MRI scan of the actual grid inside the implanted recording 
% chamber of the subject. 
%
% INPUTS:
%   GridFile:       Full path of 3D grid file (.stl/.obj/.ply/.vtk) to convert 
%   HeaderSource:   Full path of nifti (.nii) file to copy scanner
%                   alignment data from.
%   VoxelSize:      dimensions of a single voxel in mm (1x3)
%
% REQUIREMENTS (provided in ENSubfunctions folder):
%   Voxelize.m:     http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation
%   stlread.m:      http://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader/content/STLRead/stlread.m
%   Toolbox Graph:  http://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph/content/toolbox_graph/html/content.html
%   Nifti tools:    http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
%
% REVISIONS:
%   09/01/2014 - Written by Aidan Murphy (murphyap@mail.nih.gov)
%   05/02/2014 - Updated to copy scanner alignment data from existing header
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy © Copyleft 2014-2016, GNU General Public License
%==========================================================================

if nargin < 2
    GridFile     = 'Grid.stl';                                              % Specify which grid file to convert
    [file, path] = uigetfile('*.nii','Select volume to copy header from'); 	% Scan to copy scanner alignment data from
    HeaderSource = fullfile(path, file);
end
if nargin < 3
    VoxelSize   = [0.25 0.25 0.25];                                         % Set default voxel size to 0.25mm isotropic
end

addpath(genpath(fullfile(cd,'ENSubfunctions')));                    % Add required subfunctions to path
GridFormat  = GridFile(end-2:end);                                	% Get input file format
switch GridFormat                                                   % Check input format
    case 'stl'                                                      
        [v, f, n, c, stltitle] = stlread(GridFile);              	% Load .stl file
    case 'obj'
        [v,f] = read_obj(GridFile);
    case 'vtk'
        [v,f] = read_vtk(GridFile);
    case 'ply'
        [v,f] = read_ply(GridFile);
    otherwise
        fprintf(['ERROR: input grid file must be in one of the following formats:\n',...
                '\t- Stereolithography (.STL)\n'...
                '\t- Wavefront Object (.OBJ)\n'...
                '\t- Polygon File Format (.PLY)\n'...
                '\t- Visualization Toolkit (.VTK)\n']);
    	return;
end
GridSize    = max(v')-min(v');                                             	% Get grid dimensions
DimSize     = GridSize.*VoxelSize;                                          % Set volume grid size
Origin      = DimSize/2;                                                    % Find origin
Origin(3)   = 0;                                                            % Set origin of z axis to 0
gridOUTPUT  = VOXELISE(DimSize(1),DimSize(2),DimSize(3),GridFile,'xyz');    
Gridnii     = make_nii(double(gridOUTPUT),VoxelSize,Origin,[],'Grid');      % Create a new .nii file
Brainnii    = load_nii_hdr(HeaderSource);                                   % Load selected volume header
Gridnii.hdr.hist = hdr.hist;                                                % Copy scanner alignment to header 
[path, file]  = fileparts(GridFile);                                        
OutputFile  = fullfile(path, [file, '.nii']);                               % Create output file name
save_nii(Gridnii,OutputFile);                                               % Save volume as .nii