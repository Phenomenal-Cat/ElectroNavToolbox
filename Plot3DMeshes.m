function [MeshStruct] = Plot3DMeshes(MeshFiles)

%========================= Plot3DMeshes.m =================================
% 
%
%==========================================================================

if nargin == 0
    SubjectDir = 'P:/murphya/EN_data/Subjects/Layla';
    if ismac, SubjectDir = fullfile('/Volumes',SubjectDir); end
    VolumeFiles = wildcardsearch(fullfile(SubjectDir,'StructureVolumes'),'*.nii');
    MeshFiles = wildcardsearch(fullfile(SubjectDir,'VTKs'),'*.vtk');
end
   
Cmap = hsv;
facecolors = Cmap(round(linspace(1,64,numel(MeshFiles))),:);
alpha = 0.25;           


%================= LOAD MESH DATA
RevStr = '';
for i = 1:numel(MeshFiles)
    [Path, Name, temp] = fileparts(MeshFiles{i});
    MeshStruct.Names{i} = Name;
%     Msg = sprintf('%sLoading mesh %s (%d of %d)...', RevStr, MeshStruct(i).Name, i, numel(MeshFiles));
%     fprintf([RevStr, Msg]);
%   	RevStr = repmat('\b',1,length(Msg));
    [vertex,face] = read_vtk(MeshFiles{i});
    [FV(i).vertices, FV(i).faces]= read_vtk(MeshFiles{i});
    handles(i) = patch('vertices',FV(i).vertices', 'faces',FV(i).faces','edgecolor','none','facecolor',facecolors(i,:),'facealpha',alpha, 'edgealpha', alpha);
%    handles(i).PickableParts = 'none';                  % Make surface invisible to mouse clicks so that interior data points can be selected
    
    MeshStruct.Handles{i} = handles(i);
    MeshStruct.Colors(i,:) = facecolors(i,:);
    MeshStruct.Opacity(i) = alpha;
    hold on;
end
set(gca,'DataAspectRatio',[1 1 1]);
pbaspect([1 1 1]);  