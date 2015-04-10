function [fh] = ElectrodeDataToSlice(Dates, Data, Subject, SlicePlane, DataName)

%======================== ElectrodeDataToSlice.m ==========================
% This function plots a data variable onto corresponding MRI slices for the
% specified animal, for experimental session dates specified (Dates).
% There are various adjustable parameters in the top section of the
% function, but by default the function produces 3 figures:
%       1) a slice key showing the location of the slices presented in 3.
%       2) a 3D plot showing the distribution of data on atlas structures
%       3) N MRI slices in specified plane (SlicePlane) corresponding to 
%          those shown in figure 1. Each slice is overlaid with data, and
%          optionally atlas stuctures.
%
% INPUTS:
%   Dates:      a 1 x d cell array containing date strings in the format 
%               'DD-MMM-YYYY' for the data to be plotted.
%
%   Data:       a d x n matrix containing numerical values of some variable
%               (e.g. correlation coefficient) that will be plotted on the
%               color scale. Each row represents date d provided in 'Dates'
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
%   SlicePlane: 1 = corornal; 2 = axial; 3 = sagittal;
%
%   DataName:   String specifying data value type and units
%
% EXAMPLES:
%  load('LatencyData_20141016.mat')
%  ElectrodeDataToSlice(Dates, Latency, 'Layla', 3, 'RF Latency (ms)')
%
%  Dates = {'20140522','20140703','20140708','20140717','20140722','20140724','20140731'};
%  Data = repmat([1:numel(Dates)]',[1,24]);
%
%  Subject = 'Layla';
%  SlicePlane = 3;
%  ElectrodeDataToSlice(Dates, Data, 'Layla',3,'Session date')
%
% REVISIONS:
%   12/10/2014 - Written by Aidan Murphy (murphyap@mail.nih.gov)
%   15/10/2014 - Updated to transform atlas structures to native space
%   16/10/2014 - Coronal slice view added
%==========================================================================

if nargin == 0
    load('New_latency_data.mat');
    Subject = 'Layla';
	SlicePlane = 1;
    DataName = 'Latency (ms)';
end

PlotDates = 0;
if isunix && ~ismac
    opengl software;
end

%======================== ADD NECESSARY SUBFUNCTION DIRECTORIES TO PATH
Mpath = mfilename('fullpath');
Volpath = Mpath(1:strfind(Mpath,'murphy')-1);
if ~exist('load_nii','file')
    NiftiToolboxDir = fullfile(Volpath,'Toolboxes','NiftiToolbox');
    addpath(genpath(NiftiToolboxDir));
end
if ~exist('tight_subplot','file')
    APMSubDir = fullfile(Volpath,'APMSubfunctions');
	addpath(genpath(APMSubDir));
end
if ~exist('rgb','file')
	addpath(genpath(cd));
end

%========================= SET DEFAULT PATHS FOR SUBJECT
if ~exist('SlicePlane','var')
    SlicePlane = 3;       	% 1 = corornal; 2 = axial; 3 = sagittal;
end

switch Subject
    case 'Layla'
        Excel = fullfile(cd,'Subjects/Layla/LaylaElectrodeLocations.xls');              % Default Excel file
        Xform = 'Subjects/Layla/ManualXform.mat';
        MRIFile = 'Subjects/Layla/Layla_GridScan_ACPC.nii';
%         AtlasFile = 'Atlases/NeuroMaps/inia19-NeuroMaps.nii';
        AtlasFile = 'Subjects/Layla/warped_Layla_INIA19.nii';
        IndexFile = 'Subjects/Layla/warped_Layla_INIA19.mat';
%         AtlasFile = '/Subjects/Layla/winia19-NeuroMaps.nii';
        WarpMatrix = 'Subjects/Layla/inia19-t1-brain_sn.mat';
        Atlas = 1;
%         AtlasFile = 'Atlases/Frey/PaxinosPulvinar.nii';
    case 'Maliha'
        Excel = '/Subjects/Maliha/MalihaElectrodeLocations.xls';            
        Xform = '/Subjects/Maliha/Maliha_GridScan_ACPC.xform';
        MRIFile = '/Subjects/Maliha/Maliha_GridScan_ACPC.nii';
        AtlasFile = 'Atlases/inia19/inia19-NeuroMaps.nii';
    otherwise
        error('File paths for subject ''%s'' are not known! Please update .m file\n', Subject);
end


%========================= CONVERT GRID COORDINATES TO MRI COORDINATES
if PlotDates == 1
    Dates = []; 
    [GridCoords, Dates] = GetGridCoordinates(Dates, Excel, 1);      % Get grid coordinates of all contacts for requested dates
    Data = repmat([1:size(Dates,1)]',[1,24]);
%     Data(:) = 0.5;                                                  % Use grayscale
else
    [GridCoords, Dates] = GetGridCoordinates(Dates, Excel, 1);      % Get grid coordinates of all contacts for requested dates
end
GridCoords
GridCoords(:,4,:) = Data;                                           % Append input data to 4th column
GridCoords(:,5,:) = repmat([1:size(Dates,1)]',[1,size(Data,2)]);    % Assign each session an integer value
GridCoords = reshape(permute(GridCoords, [2,1,3]), [size(GridCoords,2),numel(GridCoords)/size(GridCoords,2)]);
DataLim = [min(Data(:)),max(Data(:))];                              % Find data value range

load(Xform);                                                        % Load transformation matrix
MRICoords = T*[GridCoords(1:3,:); ones(1, size(GridCoords,2))];   	% Transform grid coordinates to MRI coordinates
MRICoords([4,5],:) = GridCoords([4,5],:);                           % Add data values to 4th and 5th rows of MRIcoords


%% ======================== PLOTTING PARAMETERS ===========================

ScreenRes = get(0,'ScreenSize');            % Get screen resolution
FigPos = ScreenRes([1 2 3 3])/2;              % Set full screen height figure pos
BackgroundColor = [0.8 0.8 0.8];            % Set figure background color
PlotSliceKey = 1;
Plot3D = 1;

%================== SETTINGS FOR 3D VIEW                        
Save = 0;                                   % Automatically save figure image?
TransformAtlas = 0;                         % Transform atlas structures into native space for 3D plot?
smooth = 5;                                 % Smoothing kernel diameter (voxels)
thresh = 0.5;                               % -Inf = Otsu's method; Inf = median intensity
reduce = 0.15;                              % Reduce mesh complexity to X % of original 
ContactCmap = jet;                          % Colormap applied to Data input
IncludeEmptyChannels = 1;                   % Plot spatial locations of channels with no data?
EmptyChannelColor = [0 0 0];                % Color for empty channels
EmptyChannelAlpha = 0.4;                    % Opacity fo empty channels

HighlightIndx = [];
% HighlightColors = {[1 0 0]};%,[0 1 0]};    
% HighlightAlpha = 0.3;
% HighlightRad = 0.5;                         % Highlight radius (mm)
% HighlightDate = {'25-Mar-2014'};%,'17-Jul-2014'};
% HighlightChannel = {12};
% for h = 1:numel(HighlightDate)
%     DateIndx(h) = find(~cellfun(@isempty, strfind(Dates, HighlightDate{h})));
%     HighlightIndx(h) = DateIndx(h)+(HighlightChannel{h}-1)*numel(Dates);
% end

%==== General material lighting properties
Ambient = 0.3;                          
Diffuse = 0.5;
Specular = 0.4;
SpecExp = 6;
SpecCol = 1;

%=== Slices
SliceCmap = cool;                           % Colormap applied to slice position
SliceAlpha = 0.2;                           % Set slice alpha value
MLLims = [-16 0];                          	% Slice display limits (mm)
APLims = [-20, -8];                         
ISLims = [-8 8];

CollapseOutliers = 1;                       % Plot all data to nearest slice?
if SlicePlane==1                            % Coronal slices
%     SlicePos = -17:1:-11;               	% Set slice positions (mm)
    SlicePos = -16:2:-12;
    SubplotsW = 4;
    SubplotsH = 1;
elseif SlicePlane==2                        % Axial slices
    SlicePos = -6:0.5:6;                  	% Set slice positions (mm)
elseif SlicePlane==3                        % Sagittal
    SlicePos = -14:1:-2;  
 	SubplotsW = 5;
    SubplotsH = 2;
end
ColorIndx = round(linspace(1,64,numel(SlicePos)));
SliceColor = SliceCmap(ColorIndx,:);                % Get unique RGB color for each slice

%==== Atlas structures on slices
PlotBoundaries = 0;                                 % plot boundaries of atlas structures?
PlotArea = 1;                                       % Plot structures as solid areas?
StructFilter3D = 0;                                 % 0 = apply filter to 2D slices (faster); 1 = apply filter to 3D volume
StructGaussian = 0;                                 % Diameter of 3D Gaussian filter to apply to structures (voxels)
StructGaussSigma = 1;                               % SD of Gaussian filter
StructAlpha = 0.3;                                  % Alpha transparency of area plots
StructLineWidth = 1;                                % boundary line width
StructLineType = '-';                               % - = continuous; -- = dashed; : = dotted 
StructLabels = 0;                                   % Add text labels for structures?

%========= For atlas...
if Atlas == 1
    Structure = {'pulvinar'};%,'lateral_geniculate'};
    Bilateral = 0;      
    [StructIndex, StructNames] = GetStructureIndex(Structure, Bilateral);           % Get atlas indices for requested structure(s)
    StructIndex(5:end) = [];
    StructNames(5:end) = [];
%     StructNames([5:8,10]) = [];
%     StructIndex([5:8,10]) = [];
    LegendText = {'Inferior pulvinar','Lateral pulvinar','Medial pulvinar','Oral pulvinar'};%,'LGN'};
    colors = {'r','g','b','y','c'};
elseif Atlas == 2
    StructIndex = {3,1,2};
end
load(WarpMatrix);



%% ======================= PLOT FULL MRI SLICES ===========================

%===================== LOAD ATLAS VOLUME
AtlasNii = load_nii(AtlasFile);
AtlasNii.img = double(AtlasNii.img);
% AtlasNii.img = permute(AtlasNii.img, [2,1,3]);
VoxelSize = AtlasNii.hdr.dime.pixdim(2:4);
Origin = AtlasNii.hdr.hist.originator(1:3);
VolumeSize = size(AtlasNii.img);
if Atlas == 2
    AtlasNii.img(1:Origin(1),:,:) = 0;
end
if exist(IndexFile,'file')
    load(IndexFile);
   	load(IndexFile);
    StructIndx{1} = find(~cellfun(@isempty, strfind(StructNames,'inferior_pulvinar')));
    StructIndx{2} = find(~cellfun(@isempty, strfind(StructNames,'lateral_pulvinar')));
    StructIndx{3} = find(~cellfun(@isempty, strfind(StructNames,'medial_pulvinar')));
    StructIndx{4} = find(~cellfun(@isempty, strfind(StructNames,'oral_pulvinar')));
    StructIndx{5} = find(~cellfun(@isempty, strfind(StructNames,'lateral_geniculate')));
    StructIndx{6} = find(~cellfun(@isempty, strfind(StructNames,'medial_geniculate')));
    StructIndx{7} = find(~cellfun(@isempty, strfind(StructNames,'caudate')));
    StructIndx{8} = find(~cellfun(@isempty, strfind(StructNames,'brachium')));
    StructIndx{9} = find(~cellfun(@isempty, strfind(StructNames,'reticular')));
    StructIndx{10} = find(~cellfun(@isempty, strfind(StructNames,'lateral_ventricle')));
    StructNames = {'Pulvinar', 'LGN',  'MGN',   'Caudate', 'BSC',   'TRN',    'Ventricle'};
    Shift = [0, 3, 3];                              % Shift atlas (voxels)
    
    AtlasVol = size(AtlasNii.img);
    AtlasNii.img((end-Shift(1)):end,:,:) = [];                                      	% Shift atlas down
    AtlasNii.img = cat(1, zeros(Shift(1),AtlasVol(2),AtlasVol(3)), AtlasNii.img);   	% pad with zeros
    
    AtlasVol = size(AtlasNii.img);
    AtlasNii.img(:,:,1:Shift(3)) = [];                                              % Shift atlas down
    AtlasNii.img = cat(3, AtlasNii.img, zeros(AtlasVol(1),AtlasVol(2),Shift(3)));    % pad with zeros

    AtlasVol = size(AtlasNii.img);
    AtlasNii.img(:,(end-Shift(2):end),:) = [];                                               % Shift atlas forwards
    AtlasNii.img = cat(2, zeros(AtlasVol(1),Shift(2),AtlasVol(3)),AtlasNii.img);    % pad with zeros
    
end
StructColorNames = {'Red', 'Lime', 'blue','yellow','cyan','Orange','fuchsia', 'brown','purple','Teal'};
for i=1:numel(StructColorNames)
    StructColors(i,:) = rgb(StructColorNames{i});
end

%===================== LOAD MRI VOLUME
MRINii = load_nii(MRIFile);
MRIVoxelSize = MRINii.hdr.dime.pixdim(2:4);
MRIOrigin = MRINii.hdr.hist.originator(1:3);
MRIVolumeSize = size(MRINii.img);
Threshold = 15000;
MRINii.img(MRINii.img>Threshold) = Threshold;

if PlotSliceKey == 1
    f(1) = figure('Name','Slice key', 'Position', ScreenRes, 'Color', BackgroundColor);
    axh = tight_subplot(1, 2, 0.05, 0.05, 0.05);
%     SagSliceNum = 141;  % Centre grid hole/ 141
%     CorSliceNum = 140;
    
    SagSliceMM = -11.75;
    SagSliceNum = MRIOrigin(1)+(SagSliceMM/MRIVoxelSize(1));
    AtlasSliceNum = Origin(1)+(SagSliceMM/VoxelSize(1));
    

    %========== For sagittal view...
    SagSlice = flipud(rot90(squeeze(MRINii.img(SagSliceNum,:,:))));       	% Select sagittal slice image and reorient
    SagUpper = (MRIVolumeSize([2,3])-MRIOrigin([2,3]))*MRIVoxelSize(1);     % Get slice limits in mm relative to origin
    SagLower = -MRIOrigin([2,3])*MRIVoxelSize(1);
    axes(axh(2));
    MRh = imagesc([SagLower(1),SagUpper(1)],[SagLower(2),SagUpper(2)], SagSlice);
    hold on;
    
    %=============== PLOT ANATOMICAL STRUCTURES
    AtlasSliceImage = flipud(rot90(squeeze(AtlasNii.img(AtlasSliceNum,:,:))));
     
    for N = 1:5%numel(StructNames)    %================== For each anatomical structure...
        BinaryMask = ismember(AtlasSliceImage, StructIndx{N});
      	if StructFilter3D == 0 && StructGaussian > 0                            % If 2D slice filtering is requested...
            Fh = fspecial('gaussian', StructGaussian, StructGaussSigma);        % Create Gaussian filter
            StructAlphaLayer = imfilter(BinaryMask, Fh, 'replicate');           % Apply filter to 2D slice
            x = find(StructAlphaLayer(StructAlphaLayer>0 & StructAlphaLayer< 1));
        else
            StructAlphaLayer = BinaryMask;
        end
        
        if PlotArea==1              %================== Plot filled structure areas
            ColorLayer = repmat(zeros(size(StructAlphaLayer)), [1,1,3]);        
            for L = 1:3
                ColorLayer(:,:,L) = StructColors(N,L);
            end
            StructH(N) = imagesc([SagLower(1),SagUpper(1)],[SagLower(2),SagUpper(2)], ColorLayer);
            alpha(StructH(N),StructAlphaLayer*StructAlpha);
        end
        if PlotBoundaries == 1      %================== Plot structure boundaries
            B = bwboundaries(StructAlphaLayer);
            stat = regionprops(StructAlphaLayer,'Centroid');
            for k = 1:length(B)
                B{k}(:,2) = (B{k}(:,2)-Origin(2))*VoxelSize(2);     % Convert voxels to mm
                B{k}(:,1) = (B{k}(:,1)-Origin(3))*VoxelSize(3);
                StructLineH(N) = plot(B{k}(:,2),B{k}(:,1),[StructLineType,'k'],'color',StructColors(N,:),'linewidth',StructLineWidth);
                if StructLabels == 1 && k == 1
                    TextH(N) = text(stat(k).Centroid(1),stat(k).Centroid(2),StructNames{N},'backgroundcolor',StructColors(N,:));
                end
            end
        end
    end
    axis xy;
    set(gca,'DataAspectRatio',[1 1 1]);
    hold on;
    plot([0,0], ylim, '-w');
    plot(xlim, [0,0], '-w');
%     if SlicePlane==3
%         r = rectangle('Position',[APLims(1),ISLims(1),diff(APLims),diff(ISLims)],'LineWidth',2,'Edgecolor','c');
%     elseif SlicePlane==1
%         for s = 1:numel(SlicePos)
%             sh(s) = plot([SlicePos(s),SlicePos(s)], ISLims, '-k','Color',SliceColor(s,:),'LineWidth',2);
%         end
%     end
    colormap gray;
    freezeColors;
    FontSize = 16;
    xlabel('A-P (mm)','fontsize',FontSize);
    ylabel('I-S (mm)','fontsize',FontSize);
    title(sprintf('Sagittal slice (x = %.1fmm)', (SagSliceNum-MRIOrigin(1))*MRIVoxelSize(1)),'fontsize',FontSize);

    GridOrigin = T*[0; 0; 10; 1];
    GridTop = T*[0; 0; -20; 1];
%     plot(GridOrigin(2), GridOrigin(3), '.g');                                           % Plot grid origin
    plot([GridOrigin(2),GridTop(2)],[GridOrigin(3), GridTop(3)],'-c','linewidth',2);    % Plot center grid hole for checking orientation

    
    
    %% =================== For coronal view...
    CorSliceMM = -15;
    CorSliceNum = MRIOrigin(2)+(CorSliceMM/MRIVoxelSize(2));
    AtlasSliceNum = Origin(2)+(CorSliceMM/VoxelSize(2));
    CorSlice = flipud(rot90(squeeze(MRINii.img(:,CorSliceNum,:))));       	% Select sagittal slice image and reorient
    CorUpper = (MRIVolumeSize([1,3])-MRIOrigin([1,3]))*MRIVoxelSize(1);     % Get slice limits in mm relative to origin
    CorLower = -MRIOrigin([1,3])*MRIVoxelSize(1);
    
    
    axes(axh(1));
    imagesc([CorLower(1),CorUpper(1)],[CorLower(2),CorUpper(2)], CorSlice);
    hold on;
    
    %=============== PLOT ANATOMICAL STRUCTURES
    AtlasSliceImage = flipud(rot90(squeeze(AtlasNii.img(:,AtlasSliceNum,:))));
     
    for N = 1:5;%numel(StructNames)    %================== For each anatomical structure...
        BinaryMask = ismember(AtlasSliceImage, StructIndx{N});
      	if StructFilter3D == 0 && StructGaussian > 0                            % If 2D slice filtering is requested...
            Fh = fspecial('gaussian', StructGaussian, StructGaussSigma);        % Create Gaussian filter
            StructAlphaLayer = imfilter(BinaryMask, Fh, 'replicate');           % Apply filter to 2D slice
            x = find(StructAlphaLayer(StructAlphaLayer>0 & StructAlphaLayer< 1));
        else
            StructAlphaLayer = BinaryMask;
        end
        
        if PlotArea==1              %================== Plot filled structure areas
            ColorLayer = repmat(zeros(size(StructAlphaLayer)), [1,1,3]);        
            for L = 1:3
                ColorLayer(:,:,L) = StructColors(N,L);
            end
            StructH(N) = imagesc([CorLower(1),CorUpper(1)],[CorLower(2),CorUpper(2)], ColorLayer);
            alpha(StructH(N),StructAlphaLayer*StructAlpha);
        end
        if PlotBoundaries == 1      %================== Plot structure boundaries
            B = bwboundaries(StructAlphaLayer);
            stat = regionprops(StructAlphaLayer,'Centroid');
            for k = 1:length(B)
                B{k}(:,2) = (B{k}(:,2)-Origin(1))*VoxelSize(1);     % Convert voxels to mm
                B{k}(:,1) = (B{k}(:,1)-Origin(3))*VoxelSize(3);
                StructLineH(N) = plot(B{k}(:,2),B{k}(:,1),[StructLineType,'k'],'color',StructColors(N,:),'linewidth',StructLineWidth);
                if StructLabels == 1 && k == 1
                    TextH(N) = text(stat(k).Centroid(1),stat(k).Centroid(2),StructNames{N},'backgroundcolor',StructColors(N,:));
                end
            end
        end
    end
 
    axis xy;
    set(gca,'DataAspectRatio',[1 1 1]);
    hold on;
    plot([0,0], ylim, '-w');
    plot(xlim, [0,0], '-w');
    if SlicePlane==3
        for s = 1:numel(SlicePos)
            sh(s) = plot([SlicePos(s),SlicePos(s)], ISLims, '-k','Color',SliceColor(s,:),'LineWidth',2);
        end
    elseif SlicePlane==1
        r = rectangle('Position',[MLLims(1),ISLims(1),diff(MLLims),diff(ISLims)],'LineWidth',2,'Edgecolor',SliceColor(1,:));
    end
    colormap gray;
    freezeColors;
    xlabel('M-L (mm)','fontsize',FontSize);
    title(sprintf('Coronal slice (y = %.1fmm)', (CorSliceNum-MRIOrigin(2))*MRIVoxelSize(2)),'fontsize',FontSize);

    linkaxes(axh([1,2]),'y');
    axis equal;
    set(axh([1,2]),'ylim',[-20 40],'fontsize',16);
    set(axh(2),'ytick',[]);
    set(axh(1),'position',[0.05    0.0500    0.459    0.9000]);
end


%% =========================== PLOT 3D VIEW ===============================
AtlasNii.img = permute(AtlasNii.img,[2,1,3]); % 

if Plot3D == 1
    fh = figure('Name','3D view','Color',BackgroundColor, 'position', FigPos);
    for s = 1:numel(StructIndex)
        StructVol{s} = zeros(size(AtlasNii.img));
        StructVol{s}(ismember(AtlasNii.img, StructIndx{s})) = s;
        if (round(smooth) > 3)                                                          % smooth volume prior to edge extraction
            Vol = smooth3(StructVol{s},'gaussian',round(smooth));
        end
        FV = isosurface(StructVol{s},thresh);                                           % Create surface mesh from volume
        size(size(FV.vertices))
        FV.vertices = FV.vertices - repmat(Origin([1,2,3]),[size(FV.vertices,1),1]); 	% Translate relative to atlas origin (AC)
        FV.vertices = FV.vertices.*repmat(VoxelSize,[size(FV.vertices,1),1]);           % Scale voxels to mm
        if TransformAtlas == 0
            xyz = FV.vertices;
        elseif TransformAtlas == 1
            xyz = TransformCoordinates(FV.vertices, WarpMatrix, 0);                         % Transform atlas to native coordinates
            disp(max(xyz(:)))
            disp(min(xyz(:)))
            disp(max(FV.vertices(:)))
            disp(min(FV.vertices(:)))
        end

        %======== Plot stucture
        ph(s) = patch('vertices',xyz, 'faces',FV.faces,'facecolor',colors{s},'edgecolor','none','facealpha',StructAlpha);
        hold on;
    end
    set(gca, 'fontsize', 16);
    colormap(ContactCmap);
    legend(LegendText, 'Westoutside');


    EmptyChannelIndx = isnan(MRICoords(4,:));
    
    %============== Plot contacts
    for c = 1:numel(MRICoords(1,:))
        if ~isnan(MRICoords(4,c))
            ch(s) = DrawSphere(MRICoords(1:3,c), 0.1, 50, MRICoords(4,c), 1);
        elseif isnan(MRICoords(4,c))
            ch(s) = DrawSphere(MRICoords(1:3,c), 0.1, 50, EmptyChannelColor, EmptyChannelAlpha);
        end
    end
    if ~isempty(HighlightIndx)
        for h = 1:numel(HighlightIndx)
            ch(s+h) = DrawSphere(MRICoords(1:3,HighlightIndx(h)), HighlightRad, 50, HighlightColors{h}, HighlightAlpha);
        end
    end
    set(gca,'DataAspectRatio',[1 1 1]);
    grid on;
    lh = light('Position',[-1 1 0],'Style','infinite');
    material([Ambient Diffuse Specular SpecExp SpecCol]);
    if IncludeEmptyChannels == 0
        set(ch(EmptyChannelIndx),'visible','off');
    end
    xlabel('M-L','fontsize', 18);
    ylabel('A-P','fontsize', 18);
    zlabel('I-S','fontsize', 18);
    set(gca,'xlim',MLLims, 'ylim', APLims, 'zlim', ISLims);
    set(gca,'xtick',MLLims(1):2:MLLims(2), 'ytick', APLims(1):2:APLims(2));
    cbh = colorbar;
    set(get(cbh,'ylabel'),'string',DataName,'fontsize', 18);


    %============== Plot slices
    if SlicePlane > 0
        Xlim = xlim;
        Ylim = ylim;
        Zlim = zlim;
        if SlicePlane==1         % Coronal slices
            Ycoord = repmat(SlicePos,[4,1]);
            Xcoord = repmat(reshape([Xlim;Xlim],[4,1]),[1,numel(SlicePos)]);
            Zcoord = repmat([Zlim, Zlim([2,1])]',[1,numel(SlicePos)]);

        elseif SlicePlane==2   	% Axial slices
            Xcoord = repmat([Xlim, Xlim([2,1])]',[1,numel(SlicePos)]);
            Ycoord = repmat(reshape([Ylim;Ylim],[4,1]),[1,numel(SlicePos)]);
            Zcoord = repmat(SlicePos,[4,1]);

        elseif SlicePlane==3     % Sagittal
            Xcoord = repmat(SlicePos,[4,1]);
            Ycoord = repmat([Ylim, Ylim([2,1])]',[1,numel(SlicePos)]);
            Zcoord = repmat(reshape([Zlim;Zlim],[4,1]),[1,numel(SlicePos)]);

        end

    %     for s = 1:numel(SlicePos)
    %         Slices(s) = fill3(Xcoord(:,s), Ycoord(:,s), Zcoord(:,s), SliceColor(s,:),'facecolor',SliceColor(s,:),'facealpha',SliceAlpha,'edgecolor',SliceColor(s,:),'linewidth',2);
    %     end
    end
    % view(140,20);
    view(-135, 20);
end

% %============ Save image?
% if Save == 1
%     myaa(8);
%     export_fig(sprintf('Pulvinar3DRender_%d.png',SlicePlane),'-png','-transparent','-nocrop');
% end


%% ===================== PLOT INDIVIDUAL MR SLICES ========================
AtlasNii.img = permute(AtlasNii.img,[2,1,3]); % 

if SlicePlane == 3
    MRIndx = 1;
elseif SlicePlane == 1
    MRIndx = 2; 
end
for s = 1:numel(SlicePos)                           % For each slice position...
    ISI = diff(SlicePos([1,2]))/2;                  % Get inter-slice spacing...
    SliceIndx{s} = find((MRICoords(MRIndx,:)>= SlicePos(s)-ISI) &  (MRICoords(MRIndx,:)< SlicePos(s)+ISI));
end                                                 % Find indices for contacts that should be ploted on this slice
if CollapseOutliers == 1                            % Include data points greater than one interslice interval away on the end slices?
    SliceIndx{1} = [SliceIndx{1}, find(MRICoords(MRIndx,:)< SlicePos(1)-ISI)];
    SliceIndx{end} = [SliceIndx{end}, find(MRICoords(MRIndx,:)>= SlicePos(end)+ISI)];
end
EmptySlices = cellfun(@isempty,SliceIndx);          % Find slices containing no data
SliceIndx(EmptySlices) = [];                        % Exclude slices that contain no data!
SlicePos(EmptySlices) = [];                         
SliceColor(EmptySlices,:) = [];


% StructNames      = {'Pulvinar', 'LGN',  'MGN',   'Caudate', 'BSC',   'TRN',    'Ventricle', 'Optic tract', 'Hippocampus'};
% StructColorNames = {'Red',      'Lime', 'Orange','Cyan',    'blue','fuchsia',   'blue',      'yellow',      'Teal'};
% for i=1:numel(StructColorNames)
%     StructColors(i,:) = rgb(StructColorNames{i});
% end
% StructIndx{1} = [95,1095,96,1096,92,1092];
% StructIndx{2} = [130, 1130];
% StructIndx{3} = [114, 1114, 116, 1116];
% StructIndx{4} = [82, 1082];
% StructIndx{5} = [102, 1102];
% StructIndx{6} = [101, 1101];
% StructIndx{7} = [51, 1051, 73, 1073, 64, 1074, 201, 1201];
% StructIndx{8} = [128, 11280];
% StructIndx{9} = [448, 451, 452, 1448, 1451, 1452];






fh = figure('Name','Slice view','position', FigPos, 'Color', BackgroundColor);
axh = tight_subplot(SubplotsH, SubplotsW, 0.02, 0.02, 0.02);
for S = 1:numel(SlicePos)
    axes(axh(S));
    SliceNumber = MRIOrigin(MRIndx)+(SlicePos(S)/MRIVoxelSize(MRIndx));                  % Convert slice position (mm) to number (voxels) 
    AtlasSliceNum = Origin(MRIndx)+(SlicePos(S)/VoxelSize(MRIndx));
    
    if SlicePlane == 3
        Slice = flipud(rot90(squeeze(MRINii.img(SliceNumber,:,:))));
        AtlasSliceImage = flipud(rot90(squeeze(AtlasNii.img(AtlasSliceNum,:,:))));
    elseif SlicePlane == 1
        Slice = flipud(rot90(squeeze(MRINii.img(:,SliceNumber,:))));
        AtlasSliceImage = flipud(rot90(squeeze(AtlasNii.img(:,AtlasSliceNum,:))));
    end
    Slice = double(Slice)/max(double(Slice(:)));
    if SlicePlane == 3
        I(S) = imagesc([SagLower(1),SagUpper(1)],[SagLower(2),SagUpper(2)], repmat(Slice, [1,1,3]));
 	elseif SlicePlane == 1
        I(S) = imagesc([CorLower(1),CorUpper(1)],[CorLower(2),CorUpper(2)], repmat(Slice, [1,1,3]));
    end
    hold on;
    

    %=============== PLOT ANATOMICAL STRUCTURES
    for N = 1:numel(StructNames)    %================== For each anatomical structure...
        BinaryMask = ismember(AtlasSliceImage, StructIndx{N});
      	if StructFilter3D == 0 && StructGaussian > 0                            % If 2D slice filtering is requested...
            Fh = fspecial('gaussian', StructGaussian, StructGaussSigma);        % Create Gaussian filter
            StructAlphaLayer = imfilter(BinaryMask, Fh, 'replicate');           % Apply filter to 2D slice
            x = find(StructAlphaLayer(StructAlphaLayer>0 & StructAlphaLayer< 1));
        else
            StructAlphaLayer = BinaryMask;
        end
        
        if PlotArea==1              %================== Plot filled structure areas
            ColorLayer = repmat(zeros(size(StructAlphaLayer)), [1,1,3]);        
            for L = 1:3
                ColorLayer(:,:,L) = StructColors(N,L);
            end
            StructH(S,N) = imagesc([CorLower(1),CorUpper(1)],[CorLower(2),CorUpper(2)], ColorLayer);
            alpha(StructH(S,N),StructAlphaLayer*StructAlpha);
        end
        if PlotBoundaries == 1      %================== Plot structure boundaries
            B = bwboundaries(StructAlphaLayer);
            stat = regionprops(StructAlphaLayer,'Centroid');
            for k = 1:length(B)
                B{k}(:,2) = (B{k}(:,2)-Origin(1))*VoxelSize(1);     % Convert voxels to mm
                B{k}(:,1) = (B{k}(:,1)-Origin(3))*VoxelSize(3);
%                 B{k}(:,2) = (B{k}(:,2)-Origin(2))*VoxelSize(2);     % Convert voxels to mm
%                 B{k}(:,1) = (Origin(3)-B{k}(:,1))*VoxelSize(3);
                StructLineH(S,N) = plot(B{k}(:,2),B{k}(:,1),[StructLineType,'k'],'color',StructColors(N,:),'linewidth',StructLineWidth);
                if StructLabels == 1 && k == 1
                    TextH(S,N) = text(stat(k).Centroid(1),stat(k).Centroid(2),StructNames{N},'backgroundcolor',StructColors(N,:));
                end
            end
        end
    end
   	axis xy;
    
    %================== PLOT ELECTRODE CONTACT DATA
    hold on;
    if ~isempty(SliceIndx{S})
        for c = 1:numel(SliceIndx{S})
            if ~isnan(MRICoords(4,SliceIndx{S}(c)))
                if SlicePlane==3
                    FillCircle(MRICoords(2,SliceIndx{S}(c)), MRICoords(3,SliceIndx{S}(c)), 0.1, 20, MRICoords(4,SliceIndx{S}(c)), 1);
                elseif SlicePlane ==1
                    FillCircle(MRICoords(1,SliceIndx{S}(c)), MRICoords(3,SliceIndx{S}(c)), 0.1, 20, MRICoords(4,SliceIndx{S}(c)), 1); 
                end
            elseif isnan(MRICoords(4,SliceIndx{S}(c))) && IncludeEmptyChannels == 1
                if SlicePlane==3
                    FillCircle(MRICoords(2,SliceIndx{S}(c)), MRICoords(3,SliceIndx{S}(c)), 0.1, 20, EmptyChannelColor, EmptyChannelAlpha);
                elseif SlicePlane ==1
                    FillCircle(MRICoords(1,SliceIndx{S}(c)), MRICoords(3,SliceIndx{S}(c)), 0.1, 20, EmptyChannelColor, EmptyChannelAlpha); 
                end
            end
        end
    end
    set(gca,'DataAspectRatio',[1 1 1]);
    colormap(ContactCmap);
	set(gca,'clim',DataLim);
    if SlicePlane == 3
        set(gca,'ylim',ISLims, 'xlim', APLims,'fontsize',14);
        title(sprintf('M-L = %0.1fmm', SlicePos(S)),'fontsize',16,'backgroundcolor',SliceColor(S,:));
    elseif SlicePlane == 1
        set(gca,'ylim',ISLims, 'xlim', MLLims,'fontsize',14);
        title(sprintf('A-P = %0.1fmm', SlicePos(S)),'fontsize',16,'backgroundcolor',SliceColor(S,:));
    end
    
    if PlotDates == 1 && SlicePlane == 3
        legend(Dates{unique(MRICoords(4,SliceIndx{S}))}, 'location','EastOutside');
    end
    
end


if SlicePlane == 3
    NoPlots = SubplotsH*SubplotsW;
    if SubplotsH > 1
        set(axh([1:SubplotsW]),'xticklabel',[]);
    end
    InnerPlots = ~ismember(1:NoPlots, 1:SubplotsW:NoPlots);
    set(axh(InnerPlots),'yticklabel',[]);
end

set(axh, 'clim', DataLim);
axes(axh(S+1));
cbh = colorbar;
set(get(cbh, 'ylabel'), 'string',DataName, 'fontsize', 18);
axis off;
delete(axh(S+2:end));



end


function h = FillCircle(x,y,r,N,c,alpha)
    THETA=linspace(0,2*pi,N);
    RHO=ones(1,N)*r;
    [X,Y] = pol2cart(THETA,RHO);
    X=X+x;
    Y=Y+y;
    h=fill(X,Y,c,'EdgeColor','none','LineWidth',1);
    set(h,'FaceAlpha', alpha);
end

function h = DrawSphere(XYZ, rad, N, Color, Alpha)

	[X,Y,Z] = ellipsoid(XYZ(1),XYZ(2),XYZ(3),rad,rad,rad,N);
	
    if max(size(Color))==3
        h = surface(X,Y,Z, Z-1,'EdgeColor','none','FaceAlpha',Alpha,'FaceColor',Color);
    elseif max(size(Color))==1
        h = surface(X,Y,Z, repmat(Color,size(Z)),'EdgeColor','none','FaceAlpha',Alpha);
    end
end