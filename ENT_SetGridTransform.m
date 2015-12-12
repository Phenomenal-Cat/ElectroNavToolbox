
%====================== ENT_SetGridTransform.m ============================
% This function loads and displays an anatomical MRI with gadolinium-filled
% grid in place. The user may then manually adjust the position and
% orientation of the grid to match the scan, and the resulting
% transformation matrix is saved.
%
%
%==========================================================================

function ENT_SetGridTransform(NiiFile)

global Fig Grid nii

%============= LOAD DATA
if nargin == 0
    % [file, path]        = uigetfile({'*.nii'}, 'Select ACPC-aligned grid scan');
    % NiiFile             = fullfile(path, file);
    NiiFile             = '/Volumes/rawdata/murphya/MRI/Dexter/20151203/Dexter_20151203_MDEFT_hi_ACPC.nii';
end
nii                 = load_nii(NiiFile);                            % Load anatomical volume
nii.OriginVox     	= nii.hdr.hist.originator(1:3); 
nii.VoxDim          = nii.hdr.dime.pixdim(2:4);
nii.DimVox          = nii.hdr.dime.dim(2:4);
nii.DimMm           = nii.VoxDim.*nii.DimVox;
nii.OriginMm        = nii.OriginVox.*nii.VoxDim;
nii.AxLims          = [-nii.OriginMm; nii.DimMm-nii.OriginMm];
nii.sform           = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z];


Grid.nii              	= load_nii('Grid1.nii');                    % Load grid volume
[FV.vertices, FV.faces]	= stlread('Grid1.stl');                     % Load grid surface


Grid.rot                = [-15, 0, 0];                              % Rotations about cardinal axes (degrees)
Grid.Trans              = [7, -15, 15];                             % Translations relative to AC origin
Grid.Tform             	= makehgtform('xrotate',Grid.rot(1),'yrotate',Grid.rot(2),'zrotate',Grid.rot(3),'translate',Grid.Trans(1), Grid.Trans(2), Grid.Trans(3));
FV(2).vertices          = ENT_ApplyTform(Grid.Tform, FV(1).vertices);
Grid.Colors            	= [1 0 0; 0 1 0];


%% ============= OPEN FIGURE WINDOW
Fig.Handle      = figure('position', get(0, 'ScreenSize'),'name', 'ENT_SetGridTransform');
Fig.axh         = tight_subplot(2,3,0.05, 0.05, 0.05);
delete(Fig.axh([3,6]));
Fig.axh([3,6])  = [];









%============= PLOT 3D VIEW
axes(Fig.axh(1));
ph(1) = patch(nii.AxLims([1,1,2,2],1)', nii.AxLims([1,2,2,1],2)', zeros(1,4), zeros(1,4));%, 'Cdata', squeeze(nii.img(:,:,nii.OriginVox(3))),'FaceColor','texturemap');
ph(2) = patch(zeros(1,4), nii.AxLims([1,1,2,2],2)', nii.AxLims([1,2,2,1],3)', zeros(1,4));%, 'Cdata', squeeze(nii.img(nii.OriginVox(1),:,:)),'FaceColor','texturemap');
ph(3) = patch(nii.AxLims([1,1,2,2],1)', zeros(1,4), nii.AxLims([1,2,2,1],3)', zeros(1,4));%, 'Cdata', squeeze(nii.img(:,nii.OriginVox(2),:)),'FaceColor','texturemap');

patchTexture(ph(1), uint8(squeeze(nii.img(:,:,nii.OriginVox(3)))));
% patchTexture(ph(2), uint8(squeeze(nii.img(nii.OriginVox(1),:,:))));
% patchTexture(ph(3), uint8(squeeze(nii.img(:,nii.OriginVox(2),:))));

for g = 1:2
    Grid.h(g) = patch(FV(g), 'facecolor', Grid.Colors(g,:), 'edgecolor', 'none');     % Plot grid at origin
end
tform =	hgtransform('Parent', Fig.axh(1));
set(Grid.h(2), 'Parent', tform);
set(tform, 'Matrix', Grid.Tform);

camlight('infinite');
xlabel('X (mm)', 'fontsize', 14);
ylabel('Y (mm)', 'fontsize', 14);
zlabel('Z (mm)', 'fontsize', 14);
axis equal tight;
grid on;
view(45, 25);

%============= PLOT SLICE VIEWS
axes(Fig.axh(2));
Fig.imh(2) = imagesc(nii.AxLims(:,2)', nii.AxLims(:,1)', squeeze(nii.img(:,:,nii.OriginVox(3))));
hold on;
plot([0 0], ylim, '--w');
plot(xlim, [0 0], '--w');
xlabel('Y (mm)', 'fontsize', 14);
ylabel('X (mm)', 'fontsize', 14);
title('Axial', 'fontsize', 16);
axis equal tight;

axes(Fig.axh(3));
Fig.imh(3) = imagesc(nii.AxLims(:,3)', nii.AxLims(:,2)', squeeze(nii.img(nii.OriginVox(1),:,:)));
hold on;
plot([0 0], ylim, '--w');
plot(xlim, [0 0], '--w');
xlabel('Y (mm)', 'fontsize', 14);
ylabel('Z (mm)', 'fontsize', 14);
title('Sagittal', 'fontsize', 16);
axis equal tight;

axes(Fig.axh(4));
Fig.imh(4) = imagesc(nii.AxLims(:,1)', nii.AxLims(:,3)', squeeze(nii.img(:,nii.OriginVox(2),:)));
hold on;
plot([0 0], ylim, '--w');
plot(xlim, [0 0], '--w');
xlabel('X (mm)', 'fontsize', 14);
ylabel('Z (mm)', 'fontsize', 14);
title('Coronal', 'fontsize', 16);
axis equal tight;

colormap gray



end

%% ========================== SUBFUNCTIONS ================================
function UpdateSlices(XYZ)
global Fig Grid nii

    set(Fig.imh(2), 'cdata', squeeze(nii.img(:,:,Z)));
    set(Fig.imh(3), 'cdata', squeeze(nii.img(X,:,:)));
    set(Fig.imh(4), 'cdata', squeeze(nii.img(:,Y,:)));
    
end