function fh = RenderVolumes(VolumeDir)


if nargin==0
    SubjectID   = 'Dexter';
    Defaults    = ENT_LoadDefaults(SubjectID);
    VolumeDir   = fullfile(Defaults.ExpDir, 'Structures');
end


%=========== Set appearance defaults
smooth      = 0;
thresh      = 0.5;
alpha       = 0.25;
Ambient     = 0.3;                          
Diffuse     = 0.5;
Specular    = 0.4;
SpecExp     = 6;
SpecCol     = 1;


%=========== Load and plot structures
[filenames, pathname, filt] = uigetfile('*.nii', 'Select volumes to render', VolumeDir, 'multiselect','on');
Cmap    = jet;
Colors  = Cmap(round(linspace(1,64,numel(filenames))),:);        

for s = 1:numel(filenames)
    fprintf('Loading %s...\n', filenames{s});
    nii         = load_nii(fullfile(pathname, filenames{s}));
    Origin      = nii.hdr.hist.originator(1:3);
    VoxelSize   = nii.hdr.dime.pixdim(2:4);
    if (round(smooth) > 3)                                                          % smooth volume prior to edge extraction
        nii.img = smooth3(nii.img,'gaussian',round(smooth));
    end
    FV          = isosurface(nii.img,thresh);                                     	% Create surface mesh from volume
    FV.vertices = FV.vertices - repmat(Origin([2,1,3]),[size(FV.vertices,1),1]); 	% Translate relative to atlas origin (AC)
    FV.vertices = FV.vertices.*repmat(VoxelSize,[size(FV.vertices,1),1]);           % Scale voxels to mm
    ph(s)       = patch('vertices',FV.vertices, 'faces',FV.faces,'facecolor',Colors(s,:),'edgecolor','none','facealpha',alpha);
    hold on;
    StructNames{s} = filenames{s}(1:end-4);
    StructNames{s}(regexp(StructNames{s},'_'))= ' ';
end

%=========== Set axes appearance
set(gca,'DataAspectRatio',[1 1 1]);
grid on;
lh = light('Position',[-1 1 0],'Style','infinite');
material([Ambient Diffuse Specular SpecExp SpecCol]);
xlabel('A-P');
ylabel('M-L');
zlabel('I-S');
legend(StructNames, 'location','Eastoutside');
fh = gcf;
set(ph(1),'facealpha',1);
