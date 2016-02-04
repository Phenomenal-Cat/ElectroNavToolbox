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
    Thresh      = 15000;                        
    nii.img(nii.img>Thresh) = Thresh;
    
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

    Anatomy = permute(Anatomy,[2,1,3]);     % Switch x and y to prepare for image plotting
end


%================ CREATE STRUCTURE OUTLINES FROM MESHES
if ismember(Type, [2,3])
    StructVolumes   = wildcardsearch(Defaults.VTKdir, '*ROI.nii');
    if isempty(StructVolumes)
        error('No .nii structure volumes found in %s!', Defaults.VTKdir);
    else
        for m = 1:numel(StructVolumes)
            StructNii           = load_nii(StructVolumes{m});
            StructOriginVox     = StructNii.hdr.hist.originator(1:3);
            StructVoxSize       = StructNii.hdr.dime.pixdim(2:4);
            StructVolSize       = StructNii.hdr.dime.dim(2:4);
            StructOriginMM      = StructOriginVox.*StructVoxSize;
            for i = 1:3
                StructAxesLimsVox(i,:) = StructOriginVox(i)+(AxesLims(i,:)/StructVoxSize(i));
                StructAxRange{i} = StructAxesLimsVox(i,1):StructAxesLimsVox(i,2);
            end
%             if size(StructNii.img) ~= size(nii.img)
%                 error('MRI volume %s and structure volume %s are different sizes! (%d) Please reslice volumes to match.', Defaults.MRI, StructVolumes{m});
%             end
            for S = 1:numel(SlicePos)
                StructSliceIndx(S) = round(SlicePos(S)/StructVoxSize(Plane))+StructOriginVox(Plane);  
                switch Plane
                    case 1
                        Structures{m}(:,:,S) = squeeze(StructNii.img(StructSliceIndx(S), StructAxRange{2}, StructAxRange{3}));
                    case 2
                        Structures{m}(:,:,S) = squeeze(StructNii.img(StructAxRange{1}, StructSliceIndx(S), StructAxRange{3}));
                    case 3
                        Structures{m}(:,:,S) = squeeze(StructNii.img(StructAxRange{1}, StructAxRange{2}, StructSliceIndx(S)));
                    otherwise
                        error('Input argument ''Plane'' must be integer value in range 1-3!')
                end
            end
             Structures{m} = permute( Structures{m},[2,1,3]);     % Switch x and y to prepare for image plotting
        end
        
    end
    
   
end


%================ PLOT DATA
if TestPlot == 1
    figure;
    axh             = tight_subplot(1, numel(SlicePos), 0.02, 0.02, 0.02);
    StructColors    = jet(numel(Structures)+1);
    MaskAlpha       = 0.2; 
    FillStruct      = 0;
    
    for S = 1:numel(SlicePos)
        axes(axh(S));
        AnatomyRGB{S} = repmat(double(Anatomy(:,:,S))/double(max(max(Anatomy(:,:,S)))), [1,1,3]);
        imh(S) = image(AnatomyRGB{S});
        hold on;
        
        %================== Plot structure overlay
        if exist('Structures','var')
            for N = 1:numel(Structures)
                if FillStruct == 1
                    imsh(S,N) = imagesc(Structures{N}(:,:,S)*N,'alphadata', double(Structures{N}(:,:,S))*MaskAlpha);  	% Plot structure filled
                end
                B = bwboundaries(Structures{N}(:,:,S));                                                             % Get structure outline for each slice
                for k = 1:length(B)
%                     B{k}(:,1) = B{k}(:,1)*VoxSize(3);% + OriginMM(3);
%                     B{k}(:,2) = B{k}(:,2)*VoxSize(1);% + OriginMM(1);
                    StructLineH(S,N) = plot(B{k}(:,2), B{k}(:,1),'-w','linewidth',1);     % Plot structure boundary outline
                    set(StructLineH(S,N), 'color',StructColors(N,:));
                end
            end
        end
        
        %================== Scale images to mm
%         switch Plane
%             case 1
%                 set(imh(S),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
%                 set(imsh(S,:),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
%             case 2
%                 set(imh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
%                 set(imsh(S,:),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
%             case 3
%                 set(imh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
%                 set(imsh(S,:),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
%         end

        axis equal tight xy
        title(sprintf('%.1f mm', SlicePos(S)));
    end
    colormap jet;
end