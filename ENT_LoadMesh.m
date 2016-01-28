%========================== ENT_LoadMesh.m ================================
% This function loads 3D surface mesh files specified by the input cell array
% or string, provided they are in one of the recognized formats: .stl,
% .obj, or .vtk. The output structure FV contains vertices and faces fields
% for plotting usiong Matlab's patch.m function.
%
% EXAMPLE
%   FV = ENT_LoadMesh('Grid1.stl');
%   patch(FV,'facecolor',[1 0 0]);
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2016, GNU General Public License
%==========================================================================

function FV = ENT_LoadMesh(MeshFiles)


%============= CHECK INPUTS
if nargin == 0
    MeshFiles = uipickfiles('FilterSpec',fullfile(cd,'*.stl'));
end
if ischar(MeshFiles)
    MeshFiles = {MeshFiles};
end

%============= LOAD MESHES
for i = 1:numel(MeshFiles)
    [~,~,MeshFormat] = fileparts(MeshFiles{i});
    switch MeshFormat
        case '.vtk'
            [vertex,face] = read_vtk(MeshFiles{i});
            [FV(i).vertices, FV(i).faces]= read_vtk(MeshFiles{i});
        case '.stl'
            [vertex, face, n, c, stltitle] = stlread(MeshFiles{i});
            [FV(i).vertices, FV(i).faces]= stlread(MeshFiles{i});
        case '.obj'
            [vertex,face] = read_obj(MeshFiles{i});
            [FV(i).vertices, FV(i).faces]= read_obj(MeshFiles{i});
        otherwise
            error(sprintf('Selected mesh file format ''%s'' not recognized!',MeshFormat));
    end
    if find(size(FV(i).vertices)==max(size(FV(i).vertices)))==2
        FV(i).vertices = FV(i).vertices';
    end
    if find(size(FV(i).faces)==max(size(FV(i).faces)))==2
        FV(i).faces = FV(i).faces';
    end
end