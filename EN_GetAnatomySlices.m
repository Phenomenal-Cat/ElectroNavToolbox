function [Anatomy, Structures] = EN_GetAnatomySlices(SubjectID, Plane, SlicePos, AxesLims, Type, TestPlot)

%========================== EN_GetAnatomySlices.m =========================
% This function returns slices of anatomical data for a given subject
% within a range specified by the input parameters. The output data can be
% from either an MRI volume, outlines of anatomical structures, or both.
%
% INPUTS:
%   SubjectID:      String contaioning subject ID
%   Plane:          1 = sagittal; 2 = coronal; 3 = axial
%   SlicePos:       Vector of slice positions (mm relative to anterior comissure)
%   AxesLims:       3x2 matrix containing lower and upper limits on each row
%                   (in mm relative to the anterior commisure), for each
%                   plane.
%   Type:           1 = MRI; 2 = structure outlines; 3 = both;
%   TestPlot:       set to 1 to plot the anatomy images in a new figure.
%
% OUTPUT:
%   Anatomy:     	a 3D matrix of size XxYxN, where N is the number of
%                   slices requested and X and Y are the number of voxels
%                   within the axes limits requested. 
%   Structures:     a 1xM cell array, with each cell containing a 3D matric
%                   of size XxYxN, where M is sthe number of structure
%                   volumes, and N is the number of slices requested.
%
% EXAMPLE:
%   AxesLims = [0, 16; -30, 10; -20, 20];
%   [Anatomy, Structures] = EN_GetAnatomySlices('Dexter', 2, -18:0.5:-12, AxesLims, 3, 1);
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy © Copyleft 2014-2016, GNU General Public License
%==========================================================================

%================ CHECK INPUTS
if (size(AxesLims)~= [3,2])
    error('AxesLims input must be a 3 row by 2 column matrix')
end
Defaults    = ENT_LoadDefaults(SubjectID);


%================ GET REQUESTED SLICES FROM MRI
if ismember(Type, [1,3])
	nii         = load_nii(Defaults.MRI);
    OriginVox  	= nii.hdr.hist.originator(1:3);
    VoxSize     = nii.hdr.dime.pixdim(2:4);
    VolSize     = nii.hdr.dime.dim(2:4);
 	OriginMM    = OriginVox.*VoxSize;
    for i = 1:3
        AxesLimsVox(i,:) = OriginVox(i)+(AxesLims(i,:)/VoxSize(i));
        AxRange{i} = AxesLimsVox(i,1):AxesLimsVox(i,2);
    end

    for S = 1:numel(SlicePos)
        SliceIndx(S) = round(SlicePos(S)/VoxSize(Plane))+OriginVox(Plane);
        switch Plane
            case 1
                Anatomy(:,:,S) = squeeze(nii.img(SliceIndx(S), AxRange{2}, AxRange{3}));
            case 2
                Anatomy(:,:,S) = squeeze(nii.img(AxRange{1}, SliceIndx(S), AxRange{3}));
            case 3
                Anatomy(:,:,S) = squeeze(nii.img(AxRange{1}, AxRange{2}, SliceIndx(S)));
            otherwise
                error('Input argument ''Plane'' must be integer value in range 1-3!')
        end
    end

    Anatomy = permute(Anatomy,[2,1,3]);     % 
end


%================ CREATE STRUCTURE OUTLINES FROM MESHES
if ismember(Type, [2,3])
    StructVolumes   = wildcardsearch(Defaults.VTKdir, '*pulv.nii');
    if isempty(StructVolumes)
        error('No .nii structure volumes found in %s!', Defaults.VTKdir);
    else
        for m = 1:numel(StructVolumes)
        StructNii = load_nii(StructVolumes{m});
            if size(StructNii.img) ~= size(nii.img)
                error('MRI volume %s and structure volume %s are different sizes!', Defaults.MRI, StructVolumes{m});
            end
            for S = 1:numel(SlicePos)
                SliceIndx(S) = round(SlicePos(S)/VoxSize(Plane))+OriginVox(Plane);  
                switch Plane
                    case 1
                        Structures{m}(:,:,S) = squeeze(StructNii.img(SliceIndx(S), AxRange{2}, AxRange{3}));
                    case 2
                        Structures{m}(:,:,S) = squeeze(StructNii.img(AxRange{1}, SliceIndx(S), AxRange{3}));
                    case 3
                        Structures{m}(:,:,S) = squeeze(StructNii.img(AxRange{1}, AxRange{2}, SliceIndx(S)));
                    otherwise
                        error('Input argument ''Plane'' must be integer value in range 1-3!')
                end
            end
             Structures{m} = permute( Structures{m},[2,1,3]);     % 
        end
        
    end
    
    
%     PulvMesh        = wildcardsearch(Defaults.VTKdir, '*pulvinar_aligned.vtk');
%     if isempty(PulvMesh)
%         error('No .vtk surfaces found in %s!', Defaults.VTKdir);
%     else
%         for m = 1:numel(PulvMesh)
%             [v, f]      = read_vtk(PulvMesh{m});
%             FV.vertices = v';
%             FV.faces    = f';
%             
%             VolumeSize = size(nii.img);
%             FV.vertices = FV.vertices./VoxSize;
%             Volume = polygon2voxel(FV, VolumeSize, 'center');
%             
% %             ENT_MakeGridVolume(PulvMesh{m}, Defaults.MRI);
%             ph(m) = patch(FV);
%             
%             
            
%             %============== 
%             for S = 1:numel(SlicePos)
%                 StructNii   = load_nii(StructVolumes{m});
%                 BinaryMask  = StructNii.img;
%                 if size(nii.img) ~= size(StructNii.img)
%                     error('MRI volume %s and structure volume %s are different sizes!', Defaults.MRI, StructVolumes{m});
%                 end
%                 
%                 if StructFilter3D == 0 && StructGaussian > 0                            % If 2D slice filtering is requested...
%                     Fh = fspecial('gaussian', StructGaussian, StructGaussSigma);        % Create Gaussian filter
%                     StructAlphaLayer = imfilter(BinaryMask, Fh, 'replicate');           % Apply filter to 2D slice
%                     x = find(StructAlphaLayer(StructAlphaLayer>0 & StructAlphaLayer< 1));
%                 else
%                     StructAlphaLayer = BinaryMask;
%                 end
% 
%                 if PlotArea==1              %================== Plot filled structure areas
%                     ColorLayer = repmat(zeros(size(StructAlphaLayer)), [1,1,3]);        
%                     for L = 1:3
%                         ColorLayer(:,:,L) = StructColors(N,L);
%                     end
%                     StructH(S,N) = image(ColorLayer);
%                     alpha(StructH(S,N),StructAlphaLayer*StructAlpha);
%                 end
                if PlotBoundaries == 1      %================== Plot structure boundaries
                    B = bwboundaries(StructAlphaLayer);
                    stat = regionprops(StructAlphaLayer,'Centroid');
                    for k = 1:length(B)
                        StructLineH(S,N) = plot(B{k}(:,2),B{k}(:,1),[StructLineType,'k'],'color',StructColors(N,:),'linewidth',StructLineWidth);
                    end
                end
%             end
            
            
            
            


end




%================ PLOT DATA
if TestPlot == 1
    figure;
    axh = tight_subplot(1, numel(SlicePos), 0.02, 0.02, 0.02);
    StructColors = jet(numel(Structures)+1);
    for S = 1:numel(SlicePos)
        axes(axh(S));
        AnatomyRGB{S} = repmat(double(Anatomy(:,:,S))/double(max(max(Anatomy(:,:,S)))), [1,1,3]);
        imh(S,1) = image(AnatomyRGB{S});
        hold on;
        
        %================== Plot structure overlay
        if exist('Structures','var')
            for N = 1:numel(Structures)
                imh(S,2) = imagesc(Structures{1}(:,:,S),'alphadata', double(Structures{1}(:,:,S))*0.3);         % Plot structure filled
                B = bwboundaries(Structures{1}(:,:,S));                                                         % Get structure outline for each slice
                for k = 1:length(B)
                    StructLineH(S,N) = plot(B{k}(:,2),B{k}(:,1),'-k','color',StructColors(N,:),'linewidth',2);  % Plot outline
                end
            end
        end
        switch Plane
            case 1
                set(imh(S),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
                set(imsh(S),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
            case 2
                set(imh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
                set(imsh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
            case 3
                set(imh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
                set(imsh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
        end
        axis equal tight xy
        title(sprintf('%.1f mm', SlicePos(S)));
    end
    colormap jet
end