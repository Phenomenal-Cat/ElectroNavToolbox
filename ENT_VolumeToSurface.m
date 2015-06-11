function ENT_VolumeToSurface(Volumes, SurfFormat, OutputDir, SaveSurf)

%====================== ENT_VolumeToSurface.m =============================
% This function converts 3D volumes into 3D surface meshes which can be
% optionally displayed in a figure window before being saved. Volumes are
% filtered and thresholded before being saved as surfaces.
%
% INPUTS:
%   Volumes:        full path of 3D volume (.nii/ img) or 
%   SurfFormat:     surface file format ('vtk'/'obj'/'stl')
%   OutputDir:      directory to save surface files to
%   SaveSurf:     	1 = save surfaces without plotting ; 0 = plot without saving
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

%============ CHECK INPUTS
if ~exist('Volumes', 'var')
    Volumes = uipickfiles('REFilter', '\.nii$|\.img$', 'Prompt', 'Select volume(s) to convert');
else
    if ischar(Volumes)
        Volumes = {Volumes};
    end
end
SurfFormats = {'vtk','obj','stl'};
if ~exist('SurfFormat','var') || ~ismember(SurfFormat, SurfFormats)
	[Selection,ok] = listdlg('ListString',SurfFormats,'PromptString','Select surface output format');
    SurfFormat = SurfFormats{Selection};
end
if ~exist('SaveSurf','var')
    SaveSurf = 1;
end
if SaveSurf==1 & ~exist('OutputDir','var')
    OutputDir = uigetdir([],'Select surface output directory');
end

%=========== SET SURFACE PARAMETERS
smooth = 0;                 % Smoothing kernel
thresh = 0.5;               % Surface threshold intensity
if SaveSurf==0
    alpha = 0.25;           % Transparency (plotting only)
    Cmap = jet;             
    Colors = Cmap(round(linspace(1,64,numel(Volumes))),:);
    figure;                 % Open figure window
end
        
for s = 1:numel(Volumes)
    fprintf('Loading volume: %s (%d/%d)...\n', Volumes{s}, s, numel(Volumes));      
	SurfFilename{s} = sprintf('%s.%s', Volumes{s}(1:end-4), SurfFormat)           	% Create filename for new surface
    nii = load_nii(Volumes{s});                                                     % Load volume
    Origin = nii.hdr.hist.originator(1:3);                                          % Find volume origin
    VoxelSize = nii.hdr.dime.pixdim(2:4);                                           % Find voxel dimensions (mm)
    if (round(smooth) > 3)                                                          % smooth volume prior to edge extraction
        nii.img = smooth3(nii.img,'gaussian',round(smooth));                        
    end
    FV = isosurface(nii.img,thresh);                                                % Create surface mesh from volume
    FV.vertices = FV.vertices - repmat(Origin([2,1,3])-1,[size(FV.vertices,1),1]); 	% Translate surface relative to atlas origin (AC)
    FV.vertices = FV.vertices.*repmat(VoxelSize,[size(FV.vertices,1),1]);           % Scale voxels to mm
    
    %========= PLOT OR SAVE SURFACES
    if SaveSurf==0
        ph(s) = patch('vertices',FV.vertices,'faces',FV.faces,'facecolor',Colors(s,:),...
                      'edgecolor','none','facealpha',alpha);
        hold on;
    elseif SaveSurf==1
        switch SurfFormat
            case 'vtk'
                write_vtk(FV.vertices,FV.faces,SurfFilename{s});
            case 'obj'
                write_obj(SurfFilename{s}, FV.vertices,FV.faces);
            case 'stl'
                surf2stl(SurfFilename{s}, FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));
        end
    end
end

%========= UPDATE PLOT APPEARANCE
if SaveSurf==0
	set(gca,'DataAspectRatio',[1 1 1]);
    grid on;
    lh = light('Position',[-1 1 0],'Style','infinite');
    Ambient = 0.3;                          
    Diffuse = 0.5;
    Specular = 0.4;
    SpecExp = 6;
    SpecCol = 1;
    material([Ambient Diffuse Specular SpecExp SpecCol]);
end