function MeshStruct = ENT_Plot3DMeshes(MeshFiles)

%========================= ENT_Plot3DMeshes.m =============================
% Plots all 3D surface meshes listed in the input array of full path
% strings 'MeshFiles'. If 'MeshFiles' is a directory, all mesh files of a
% supported format in that directory will be plotted. Supported file types
% are: .vtk; .obj; .stl; 
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, � Copyleft 2015, GNU General Public License
%==========================================================================

% addpath(genpath(cd))
if nargin == 0 || numel(MeshFiles) == 0
    MeshFiles = uipickfiles('FilterSpec',fullfile(cd,'*.stl'));
end
if ~iscell(MeshFiles) && isdir(MeshFiles)    
    MeshFiles = wildcardsearch(MeshFiles,'*.vtk');
end

Cmap            = jet;
facecolors      = Cmap(round(linspace(1,64,numel(MeshFiles))),:);
DefaultAlpha    = 0.25;       
 

%================= LOAD MESH DATA
for i = 1:numel(MeshFiles)
    [Path, Name, temp] = fileparts(MeshFiles{i});
    MeshStruct.Names{i} = Name;
	fprintf('Loading mesh %s (%d of %d)...\n', MeshStruct.Names{i}, i, numel(MeshFiles));
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
    handles(i) = patch('vertices',FV(i).vertices, 'faces',FV(i).faces,'edgecolor','none','facecolor',facecolors(i,:),'facealpha',DefaultAlpha, 'edgealpha', DefaultAlpha);
    try
        handles(i).PickableParts = 'none';                  % Make surface invisible to mouse clicks so that interior data points can be selected
    end
    MeshStruct.Handles{i}   = handles(i);
    MeshStruct.Colors(i,:)  = facecolors(i,:);
    MeshStruct.Opacity(i)   = DefaultAlpha;
    hold on;
end
set(gca,'DataAspectRatio',[1 1 1]);                     
pbaspect([1 1 1]);                
grid on;