
%============================== Plot Slice View ===========================
% This script plots anatomical structures from the specified atlas onto the
% MRI slices specified. The atlas can be in template space (e.g. the
% NeuroMaps atlas plotted onto the INIA19 MRI template) or native space
% (e.g. a spatially normalized version of the atlas on an individual's
% MRI).
%
%
%==========================================================================

UseAtlasSpace = 1;

[ENTroot, temp, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(ENTroot));
if UseAtlasSpace == 1
    Atlasfile = fullfile(ENTroot, 'Atlases/inia19/inia19-NeuroMaps.nii');
    MRIfile = fullfile(ENTroot,'Atlases/inia19/inia19-t1-brain.nii');
%   	Atlasfile = fullfile(ENTroot, 'Atlases/NeuroMaps/inia19-NeuroMaps.nii');
%     MRIfile = fullfile(ENTroot,'Atlases/NeuroMaps/inia19-t1-brain.nii');
    MRIThresholdIntensity = 180;
else
   	Atlasfile = fullfile(ENTroot, 'Subjects/Layla/warped_Layla_INIA19.nii');
   	IndexFile = fullfile(ENTroot, 'Subjects/Layla/warped_Layla_INIA19.mat');
    MRIfile = fullfile(ENTroot,'Subjects/Layla/Layla_GridScan_ACPC.nii');
    MRIThresholdIntensity = 15000;
    Shift = [2, 4, 4];                              % Shift atlas (voxels)
end

AtlasNii = load_nii(Atlasfile);                     % Load spatially normalized atlas volume
MRINii = load_nii(MRIfile);                         % Load anatomical ACPC aligned volume
AtlasNii.img = permute(AtlasNii.img,[1,3,2]);       % Permute volume to display in correct orientation
MRINii.img = permute(MRINii.img,[1,3,2]);           
AtlasVolmm = prod(size(AtlasNii.img).*AtlasNii.hdr.dime.pixdim(2:4));
MRIVolmm = prod(size(MRINii.img).*MRINii.hdr.dime.pixdim(2:4));

VM = '';
if AtlasNii.hdr.dime.pixdim(2:4) ~= MRINii.hdr.dime.pixdim(2:4)
    VM = sprintf('%s\t- voxel sizes: [%d, %d, %d] VS [%d, %d, %d];\n', VM, AtlasNii.hdr.dime.pixdim(2:4), MRINii.hdr.dime.pixdim(2:4));
end
if AtlasNii.hdr.hist.originator(1:3) ~= MRINii.hdr.hist.originator(1:3)
    VM = sprintf('%s\t- origin:      [%d, %d, %d] VS [%d, %d, %d];\n', VM, AtlasNii.hdr.hist.originator(1:3), MRINii.hdr.hist.originator(1:3));
end
if ~isempty(VM) 
     error(['The NATIVE ATLAS and NATIVE T1 volumes provided do not have the same:\n%s',...
            'Reslice one volume to match before trying again!\n'], VM);
end
% if AtlasVolmm > MRIVolmm                                % If atlas volume is larger than T1 volume
%     MRIorigin = MRINii.hdr.hist.originator(1:3);
%     MRIdistfromorigin = (size(MRINii.img)-MRIorigin);
%     Atlasorigin = AtlasNii.hdr.hist.originator(1:3);
%     Indx([1,3,5]) = abs(Atlasorigin-MRIorigin);
%     Indx([2,4,6]) = size(AtlasNii.img)-MRIdistfromorigin;
%     
%     AtlasNii.img = AtlasNii.img(Indx(1):Indx(2), Indx(3):Indx(4), Indx(5):Indx(6));
% end
     

%================== Atlas offset adjustments
if exist('Shift','var')
    AtlasVol = size(AtlasNii.img);
    AtlasNii.img((end-Shift(1)):end,:,:) = [];                                      	% Shift atlas down
    AtlasNii.img = cat(1, zeros(Shift(1),AtlasVol(2),AtlasVol(3)), AtlasNii.img);   	% pad with zeros
    
    AtlasVol = size(AtlasNii.img);
    AtlasNii.img(:,:,(end-Shift(3)-1):end) = [];                                     % Shift atlas forwards
    AtlasNii.img = cat(3, zeros(AtlasVol(1),AtlasVol(2),Shift(3)), AtlasNii.img);    % pad with zeros

    AtlasVol = size(AtlasNii.img);
    AtlasNii.img(:,1:Shift(2),:) = [];                                               % Shift atlas down
    AtlasNii.img = cat(2, AtlasNii.img, zeros(AtlasVol(1),Shift(2),AtlasVol(3)));    % pad with zeros
end


%====================== Set plot parameters
PlotBoundaries = 1;                                 % 0 = plot solid areas; 1 = plot boundaries
PlotArea = 1;
StructFilter3D = 0;                                 % 0 = apply filter to 2D slices (faster); 1 = apply filter to 3D volume
StructGaussian = 0;                                 % Diameter of 3D Gaussian filter to apply to structures (voxels)
StructGaussSigma = 1;                               % SD of Gaussian filter
StructAlpha = 0.2;                                  % Alpha transparency of area plots
StructLineWidth = 1;                                % boundary line width
StructLineType = '-';                               % - = continuous; -- = dashed; : = dotted 
StructLabels = 0;                                   % Add text labels for structures?
AxisView = 2;                                       % 1 = sag; 2 = cor; 3 = ax
if AxisView == 1
    SlicePos = -14:1:-4;                                % Set slice positions (mm)
elseif AxisView == 2
    SlicePos = -18:-10;
    AtlasNii.img = permute(AtlasNii.img,[3,2,1]);
    MRINii.img = permute(MRINii.img,[3,2,1]);
end

if UseAtlasSpace == 0
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
    StructNames      = {'inf_Pulvinar', 'l_Pulvinar','m_Pulvinar','o_Pulvinar','LGN',  'MGN',   'Caudate', 'BSC',   'TRN',    'Ventricle'};
else
    StructNames      = {'Pulvinar', 'LGN',  'MGN',   'Caudate', 'BSC',   'TRN',    'Ventricle', 'Optic tract', 'Hippocampus'};
    StructIndx{1} = [95,1095,96,1096,92,1092];
    StructIndx{2} = [130, 1130];
    StructIndx{3} = [114, 1114, 116, 1116];
    StructIndx{4} = [82, 1082];
    StructIndx{5} = [102, 1102];
    StructIndx{6} = [101, 1101];
    StructIndx{7} = [51, 1051, 73, 1073, 64, 1074, 201, 1201];
    StructIndx{8} = [128, 11280];
    StructIndx{9} = [448, 451, 452, 1448, 1451, 1452];
    % [122, 1122, 123, 1123];   % Habenula
end

PlotAmygdala = 1;
if PlotAmygdala == 1
    StructNames      = {'accessory basal','anterior amygdalar area','basal nucleus','central nucleus','cortical nucleus','lateral nucleus','medial nucleus','paralaminar nucleus','substantia_innominata'};
    StructIndx{1} = [191, 1191];
    StructIndx{2} = [496, 1496];
    StructIndx{3} = [192        1192];
    StructIndx{4} = [482        1482];
    StructIndx{5} = [498        1498];
    StructIndx{6} = [189        1189];
    StructIndx{7} = [484        1484];
    StructIndx{8} = [447        1447];
    StructIndx{9} = [292        1292];
%     StructIndx{10} = [94         188        1094        1188];
	if AxisView == 1
        SlicePos = -15:1:-5;
   	elseif AxisView == 2
    	SlicePos = -5:3;
    end
end


StructColorNames = {'Red','Lime', 'blue','yellow','Cyan','Orange','fuchsia', 'purple', 'brown',      'Teal'};
for i=1:numel(StructColorNames)
    StructColors(i,:) = rgb(StructColorNames{i});
end

MRINii.img(MRINii.img>MRIThresholdIntensity) = MRIThresholdIntensity;   % Threshold voxel intensities in MR image


%% =========================== PLOT SLICE DATA ============================
fh = figure;
Axh = tight_subplot(3, 4, 0.02, 0.02, 0.02);
for S = 1:numel(SlicePos)           %================= For each slice requested...
    axes(Axh(S));
    
    MRSliceNum(S) = SlicePos(S)/MRINii.hdr.dime.pixdim(2)+MRINii.hdr.hist.originator(AxisView);
    MRISliceImage = squeeze(MRINii.img(MRSliceNum(S),:,:));
    
    %========= Adjust contrast
    MRISliceImage(1,1) = 0; MRISliceImage(1,2) = MRIThresholdIntensity;
    
    MRIH = imagesc(MRISliceImage);                                                      
    hold on;
    
    AtlasSliceNum(S) = SlicePos(S)/AtlasNii.hdr.dime.pixdim(2)+AtlasNii.hdr.hist.originator(AxisView);
    AtlasSliceImage = squeeze(AtlasNii.img(AtlasSliceNum(S),:,:));


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
            StructH(S,N) = image(ColorLayer);
            alpha(StructH(S,N),StructAlphaLayer*StructAlpha);
        end
        if PlotBoundaries == 1      %================== Plot structure boundaries
            B = bwboundaries(StructAlphaLayer);
            stat = regionprops(StructAlphaLayer,'Centroid');
            for k = 1:length(B)
                StructLineH(S,N) = plot(B{k}(:,2),B{k}(:,1),[StructLineType,'k'],'color',StructColors(N,:),'linewidth',StructLineWidth);
                if StructLabels == 1 && k == 1
                    TextH(S,N) = text(stat(k).Centroid(1),stat(k).Centroid(2),StructNames{N},'backgroundcolor',StructColors(N,:));
                end
            end
        end
    end
    
    colormap gray;
    axis xy off;
    set(gca,'DataAspectRatio',[1 1 1],'TickDir','out');
    box off;
    title(sprintf('X = %0.1f mm', SlicePos(S)));
end


%=================== Add figure legend in spare axes
axes(Axh(S+1));
for N = 1:numel(StructNames)
    ph(N) = rectangle('position',[10,N*10,5,5],'FaceColor', StructColors(N,:));
    text(20, N*10+2, StructNames{N});
end
set(gca,'xlim',[0 100],'DataAspectRatio',[1 1 1]);
axis off;

delete(Axh(numel(SlicePos)+2:end));
Axh(numel(SlicePos)+1:end) = [];
linkaxes(Axh);
if UseAtlasSpace == 1
	set(Axh,'xlim',[70,110],'ylim',[40, 80]);       
else
    set(Axh,'ylim',[80 250],'xlim',[50 310]);
end

Origin = AtlasNii.hdr.hist.originator(1:3);
ISlim = [-10 10];
ISlimVox = (ISlim/0.25)+Origin(3);
MLlim = [-18 0];
MLlimVox = (MLlim/0.25)+Origin(1);

set(Axh,'ylim',ISlimVox,'xlim',MLlimVox);
