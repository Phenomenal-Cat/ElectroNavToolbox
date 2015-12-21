% function [out] = ElectrodeDataToSlice(Dates, Data, Subject)

%======================== ElectrodeDataToSlice.m ==========================
% 
%
%
% INPUTS:
%
%   Dates:      a 1 x d cell array containing date strings in the format 
%               'DD-MMM-YYYY' for the data to be plotted.
%
%   Data:       a d x n matrix containing numerical values of some variable
%               (e.g. correlation coefficient) that will be plotted on the
%               color scale. Each row represents a date provided in 'Dates'
%               and each column represents an electrode contact number
%               (arranged in ascending physical order, i.e. column 1 = 
%               contact closest to the electrode tip).
%
%   Subject:    String identifying subject. This determines the default
%               paths for the following:
%
%               Excel:      an Excel file containing electrode position data for all
%                           session dates.
%               Xform:      a 4 x 4 transformation matrix for converting grid-centred
%                           coordinates into native ACPC aligned MRI-space coordinate 
%                           frame.
%               MRIFile:    Native MRI grid scan of the individual subject, ACPC aligned
%                           with origin set to the anterior commissure.
%
% EXAMPLE:
%
%   Dates = {'31-Jul-2014','10-Jul-2014','04-Jun-2014'};
%   Data = rand(3,24);;
%   ElectrodeDataToSlice(Dates, Data, 'Layla')
%
% REVISIONS:
%
%   12/10/2014 - Written by Aidan Murphy (murphyap@mail.nih.gov)
%
%==========================================================================

Subject = 'Layla';

%========================= SET DEFAULT PATHS FOR SUBJECT
switch Subject
    case 'Layla'
        Excel = '/Subjects/Layla/LaylaElectrodeLocations.xls';              % Default Excel file
        Xform = '/Subjects/Layla/Layla_GridScan_ACPC.xform';
        MRIFile = '/Subjects/Layla/Layla_GridScan_ACPC.nii';
        AtlasFile = 'Atlases/NeuroMaps/inia19-NeuroMaps.nii';
        Atlas = 1;
%         AtlasFile = 'Atlases/Frey/PaxinosPulvinar.nii';
    case 'Maliha'
        Excel = '/Subjects/Maliha/MalihaElectrodeLocations.xls';            
        Xform = '/Subjects/Maliha/Maliha_GridScan_ACPC.xform';
        MRIFile = '/Subjects/Maliha/Maliha_GridScan_ACPC.nii';
        AtlasFile = 'Atlases/NeuroMaps/inia19-NeuroMaps.nii';
    otherwise
        error('File paths for subject ''%s'' are not known! Please update .m file\n', Subject);
end


%========================= CONVERT GRID COORDINATES TO MRI COORDINATES
% GridCoords = GetGridCoordinates(Dates, Excel, 1);           % Get grid coordinates of all contacts for requested dates
% T = load(Xform);                                            % Load transformation matrix
% GridCoords = reshape(permute(GridCoords, [2,1,3]), [3,numel(GridCoords)/3]);    
% MRICoords = GridCoords.*T;                                  % Transform grid coordinates to MRI coordinates


%========================= SETTINGS FOR 3D VIEW ===========================                        
Save = 0;                                   % Automatically save figure image?
smooth = 5;                                 % Smoothing kernel diameter (voxels)
thresh = 0.5;                               % -Inf = Otsu's method; Inf = median intensity
reduce = 0.15;                              % Reduce mesh complexity to X % of original 
colors = {'r','g','b','y'};
alpha = 1;

%==== General material lighting properties
Ambient = 0.3;                          
Diffuse = 0.5;
Specular = 0.4;
SpecExp = 6;
SpecCol = 1;

%=== Slices
SlicePlane = 3;                             % 0 = none; 1 = corornal; 2 = axial; 3 = sagittal;
SliceColor = 'b';                           % Set slice edge color
SliceAlpha = 0.2;                           % Set slice alpha value
MLLims = [0 16];                            % Slice display limits
APLims = [-20, -8];
ISLims = [-8 8];
if SlicePlane==1                            % Coronal slices
    SlicePos = -11:0.5:-17;               	% Set slice positions (mm)
elseif SlicePlane==2                        % Axial slices
    SlicePos = -6:0.5:6;                  	% Set slice positions (mm)
elseif SlicePlane==3                        % Sagittal
    SlicePos = 2:1:14;  
end

Atlas = 1;
if Atlas == 1
    Structure = 'pulvinar';
    Bilateral = 0;      
    [StructIndex, StructNames] = ENT_GetStructureIndex(Structure, Bilateral);           % Get atlas indices for requested structure(s)
    StructIndex(1:4) = [];
    StructNames(1:4) = [];
elseif Atlas == 2
    StructIndex = {3,1,2};
end


%================= ASSIGN CONTACT COORDINATES TO NEAREST SLICE ============
% if SlicePlane==3                        % Sagittal
%     MLcontactpos = round(MRICoords(:,SlicePos
% end



%% ======================= Plot MR slices =================================
MRINii = load_nii(MRIFile);
MRIVoxelSize = MRINii.hdr.dime.pixdim(2:4);
MRIOrigin = MRINii.hdr.hist.originator(1:3);
MRIVolumeSize = size(MRINii.img);

figure;
axh = tight_subplot(1, 3, 0.02, 0.02, 0.02);
SagSliceNum = 141;
CorSliceNum = 140;

%========== For sagittal view...
SagSlice = flipud(rot90(squeeze(MRINii.img(141,:,:))));                 % Select sagittal slice image and reorient
SagUpper = (MRIVolumeSize([2,3])-MRIOrigin([2,3]))*MRIVoxelSize(1);     % Get slice limits in mm relative to origin
SagLower = -MRIOrigin([2,3])*MRIVoxelSize(1);
axes(axh(1));
imagesc([SagLower(1),SagUpper(1)],[SagLower(2),SagUpper(2)], SagSlice);
axis xy;
set(gca,'DataAspectRatio',[1 1 1]);
hold on;
plot([0,0], ylim, '--r');
plot(xlim, [0,0], '--r');
if SlicePlane==3
    r = rectangle('Position',[APLims(1),ISLims(1),diff(APLims),diff(ISLims)],'LineWidth',2,'Edgecolor',SliceColor(1,:));
elseif SlicePlane==1
    
end
colormap gray;
xlabel('A-P (mm)');
ylabel('I-S (mm)');
title(sprintf('Sagittal slice (x = %.1fmm)', (SagSliceNum-MRIOrigin(1))*MRIVoxelSize(1)));


%========== For coronal view...
CorSlice = flipud(rot90(squeeze(MRINii.img(:,CorSliceNum,:))));       	% Select sagittal slice image and reorient
CorUpper = (MRIVolumeSize([1,3])-MRIOrigin([1,3]))*MRIVoxelSize(1);     % Get slice limits in mm relative to origin
CorLower = -MRIOrigin([1,3])*MRIVoxelSize(1);
axes(axh(2));
imagesc([CorLower(1),CorUpper(1)],[CorLower(2),CorUpper(2)], CorSlice);
axis xy;
set(gca,'ylim',[-20 50]);
set(gca,'DataAspectRatio',[1 1 1]);
hold on;
plot([0,0], ylim, '--r');
plot(xlim, [0,0], '--r');
if SlicePlane==3
    for s = 1:numel(SlicePos)
        sh(s) = plot([SlicePos(s),SlicePos(s)], ISLims, '-k','Color',SliceColor(s,:));
    end
elseif SlicePlane==1
    r = rectangle('Position',[MLLims(1),ISLims(1),diff(MLLims),diff(ISLims)],'LineWidth',2,'Edgecolor',SliceColor(1,:));
end
colormap gray;
xlabel('M-L (mm)');
title(sprintf('Coronal slice (y = %.1fmm)', (CorSliceNum-MRIOrigin(2))*MRIVoxelSize(2)));


%========== For 3D view...

%===================== LOAD VOLUMES
AtlasNii = load_nii(AtlasFile);
VoxelSize = AtlasNii.hdr.dime.pixdim(2:4);
Origin = AtlasNii.hdr.hist.originator(1:3);
VolumeSize = size(AtlasNii.img);
if Atlas == 2
    AtlasNii.img(1:Origin(1),:,:) = 0;
end
axes(axh(3));
for s = 1:numel(StructIndex)
    StructVol{s} = zeros(size(AtlasNii.img));
    StructVol{s}(ismember(AtlasNii.img, StructIndex{s})) = s;
    if (round(smooth) > 3)                                                          % smooth volume prior to edge extraction
        Vol = smooth3(StructVol{s},'gaussian',round(smooth));
    end
    FV = isosurface(StructVol{s},thresh);                                           % Create surface mesh from volume
    FV.vertices = FV.vertices - repmat(Origin([2,1,3]),[size(FV.vertices,1),1]); 	% Translate relative to atlas origin (AC)
    FV.vertices = FV.vertices.*repmat(VoxelSize,[size(FV.vertices,1),1]);           % Scale voxels to mm
    ph(s) = patch('vertices',FV.vertices, 'faces',FV.faces,'facecolor',colors{s},'edgecolor','none','facealpha',alpha);
    hold on;
end
set(gca,'DataAspectRatio',[1 1 1]);
grid on;
lh = light('Position',[-1 1 0],'Style','infinite');
material([Ambient Diffuse Specular SpecExp SpecCol]);
xlabel('A-P');
ylabel('M-L');
zlabel('I-S');
set(gca,'ylim',MLLims, 'xlim', APLims, 'zlim', ISLims);
LegendText = {'Inferior pulvinar','Lateral pulvinar','Medial pulvinar','Oral pulvinar'};
legend(LegendText, 'Eastoutside');


%============== Plot slices
if SlicePlane > 0
    Xlim = xlim;
    Ylim = ylim;
    Zlim = zlim;
    if SlicePlane==1         % Coronal slices
      	Xcoord = repmat(SlicePos,[4,1]);
        Ycoord = repmat(reshape([Ylim;Ylim],[4,1]),[1,numel(SlicePos)]);
        Zcoord = repmat([Zlim, Zlim([2,1])]',[1,numel(SlicePos)]);
        
    elseif SlicePlane==2   	% Axial slices
        Xcoord = repmat([Xlim, Xlim([2,1])]',[1,numel(SlicePos)]);
        Ycoord = repmat(reshape([Ylim;Ylim],[4,1]),[1,numel(SlicePos)]);
        Zcoord = repmat(SlicePos,[4,1]);
        
    elseif SlicePlane==3     % Sagittal
       	Xcoord = repmat([Xlim, Xlim([2,1])]',[1,numel(SlicePos)]);
        Ycoord = repmat(SlicePos,[4,1]);
        Zcoord = repmat(reshape([Zlim;Zlim],[4,1]),[1,numel(SlicePos)]);
        
    end
    
    cmap = jet;
    ColorIndx = round(linspace(1,64,numel(SlicePos)));
    SliceColor = cmap(ColorIndx,:);
    for s = 1:numel(SlicePos)
        Slices(s) = fill3(Xcoord(:,s), Ycoord(:,s), Zcoord(:,s), SliceColor(s,:),'facecolor',SliceColor(s,:),'facealpha',SliceAlpha,'edgecolor',SliceColor(s,:),'linewidth',2);
    end
end
% view(140,20);
view(-135, 20);


linkaxes(axh([1,2]),'y');
set(axh([1,2]),'ylim',[-20 50]);
axis equal;



%============ Save image?
if Save == 1
    myaa(8);
    export_fig(sprintf('Pulvinar3DRender_%d.png',SlicePlane),'-png','-transparent','-nocrop');
end




%% ===================== PLOT INDIVIDUAL MR SLICES ========================
figure;
axh = tight_subplot(3, 5, 0.02, 0.02, 0.02);
for S = 1:numel(SlicePos)
    axes(axh(S+1));
    SliceNumber = (-SlicePos(S)/MRIVoxelSize(1))+MRIOrigin(1);   % Convert slice position (mm) to number (voxels) 
    Slice = squeeze(MRINii.img(SliceNumber,:,:));
    I = imagesc(rot90(Slice));
    axis xy;
    set(gca,'DataAspectRatio',[1 1 1]);
    colormap gray
    
    set(gca,'ylim',ISLims, 'xlim', APLims);
    title(sprintf('M-L = %0.1fmm', -SlicePos(S)));
    axis off;
end
delete(axh(S+2:end));
