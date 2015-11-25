
%========================== ENT_DataCloudToSlice.m ========================
% This function takes data points within a 3D volume and assigns them to
% slice planes. For each slice
%
% INPUTS:   CloudData.XYZ:      n x 3 matrix containing 3D coordinates for 
%                               each data point (mm).
%           CloudData.Values:   data value for each data point.
%           CloudData.Weights:  optional vector of weightings (0-1) for
%                               each point.
%
%           SliceParams.Plane:  1 = sagital; 2 = coronal; 3 = axial
%           SliceParams.ISD:    inter-slice distance (mm)
%           SliceParams.Interp: 1 = nearest neighbour; 2 = trilinear;
%           SliceParams.Filt:   1 = 
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function SliceData = ENT_DataCloudToSlice(CloudData, SliceParams)
global SliceParams InPlane

if nargin==0
    
end

DataFile = '/Volumes/PROJECTS/murphya/Physio/MapPlotTools/MapData/Dexter_PeakToTroughDuration_S=20_N=1215';
load(DataFile)
CloudData       = Contact;
CloudData.XYZ2   = squeeze(reshape(CloudData.XYZ, [size(CloudData.XYZ,1)*size(CloudData.XYZ,2), 1, 3]));


SliceParams.Plane       = 2;
SliceParams.ISD         = 2;
SliceParams.GridRes     = 0.25;          
SliceParams.Interp      = 1;
SliceParams.Filt        = 0;
SliceParams.Range       = [floor(min(CloudData.XYZ2)); ceil(max(CloudData.XYZ2))];
SliceParams.SliceEdges  = SliceParams.Range(1,SliceParams.Plane):SliceParams.ISD:SliceParams.Range(2,SliceParams.Plane);
SliceParams.SlicePos    = SliceParams.SliceEdges(1:end-1)+(SliceParams.ISD/2);
SliceParams.Filt        = fspecial('gaussian', [5, 5], 1.5);
SliceParams.Filt2     	= fspecial('gaussian', [3, 3], 1);

PlaneIDs                = {'X','Y','Z'};
PlaneLabels           	= {'Medial-lateral (mm)','Posterior-anterior (mm)','Inferior-superior (mm)'};
InPlane                 = [1,2,3];
InPlane(SliceParams.Plane) = [];

figure;
axh = tight_subplot(2, numel(SliceParams.SlicePos), 0.05, 0.05, 0.05);

for s = 1:numel(SliceParams.SlicePos)
    PosIndices{s}   = find(CloudData.XYZ2(:,SliceParams.Plane)>= SliceParams.SliceEdges(s) & CloudData.XYZ2(:,SliceParams.Plane)< SliceParams.SliceEdges(s+1));
    ChannelIndx{s} 	= ceil(PosIndices{s}/size(CloudData.XYZ,1));
    DateIndx{s}    	= PosIndices{s}-(ChannelIndx{s}-1)*size(CloudData.XYZ,1);
    DataIndices{s}  = find(ismember(CloudData.CellIndxData(:,[2,3]), [DateIndx{s}, ChannelIndx{s}],'rows'));
    
    for ch = 1:numel(ChannelIndx{s})
        ChIndx = find(ismember(CloudData.CellIndxData(:,[2,3]), [DateIndx{s}(ch), ChannelIndx{s}(ch)],'rows'));
        ChannelAverages{s}(ch) = mean(CloudData.ColorVals(ChIndx)); 
        ChannelStd{s}(ch) = std(CloudData.ColorVals(ChIndx));
    end
    
    %============ Plot raw contact positions
    axes(axh(s));
    plot(CloudData.XYZ2(PosIndices{s}, InPlane(1)), CloudData.XYZ2(PosIndices{s}, InPlane(2)),'.k');
    title(sprintf('%s = %.1f mm', PlaneIDs{SliceParams.Plane}, SliceParams.SlicePos(s)),'fontsize',18);
    xlabel(PlaneLabels{InPlane(1)},'fontsize',15);
    ylabel(PlaneLabels{InPlane(2)},'fontsize',15);
    grid on
    
    %============ Plot interpolated data
    axes(axh(s)+numel(SliceParams.SlicePos));
    Grid(s).Data = zeros(diff(SliceParams.Range(:,InPlane(2)))/SliceParams.GridRes, diff(SliceParams.Range(:,InPlane(1)))/SliceParams.GridRes);
    Grid(s).Alpha = Grid(s).Data;
    X = floor((CloudData.XYZ2(PosIndices{s}, InPlane(1)) - SliceParams.Range(1,InPlane(1)))/SliceParams.GridRes);
    Y = floor((CloudData.XYZ2(PosIndices{s}, InPlane(2)) - SliceParams.Range(1,InPlane(2)))/SliceParams.GridRes);
    for p = 1:numel(X)
        Grid(s).Data(Y(p), X(p)) = ChannelAverages{s}(p);
        Grid(s).Alpha(Y(p), X(p)) = 1;
    end
    Grid(s).DataFilt = imfilter(Grid(s).Data, SliceParams.Filt, 0);
    Grid(s).AlphaFilt = imfilter(Grid(s).Alpha, SliceParams.Filt2, 0);
    Imh(s) = imagesc(SliceParams.Range(1,InPlane(1)):SliceParams.Range(2,InPlane(1)), SliceParams.Range(1,InPlane(2)):SliceParams.Range(2,InPlane(2)), Grid(s).DataFilt);
    set(Imh(s), 'alphadata', Grid(s).AlphaFilt);
    axis xy
    xlabel(PlaneLabels{InPlane(1)},'fontsize',15);
    ylabel(PlaneLabels{InPlane(2)},'fontsize',15);
    grid on
    
end
linkaxes(axh);
set(axh, 'xlim', SliceParams.Range(:,InPlane(1)), 'ylim', SliceParams.Range(:,InPlane(2)), 'color', [0.5 0.5 0.5]);
colormap jet
colorbar

end
