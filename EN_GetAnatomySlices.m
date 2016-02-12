function [Anatomy, Structures, Outlines] = EN_GetAnatomySlices(SubjectID, Plane, SlicePos, AxesLims, Type, TestPlot)

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
%   Structures:     a 1xM cell array, with each cell containing a 3D matrix
%                   of size XxYxN, where M is the number of structure
%                   volumes, and N is the number of slices requested.
%   Outlines:       an SxMxk cell array, where M is the number of
%                   tructures, S is the number of slices and k is the number
%                   of separate parts of the structure. Each cell contains a 2xN
%                   matrix of in-slice coordinates.
%
% EXAMPLE:
%   AxesLims = [0, 16; -30, 10; -20, 20];
%   [Anatomy, Structures] = EN_GetAnatomySlices('Dexter', 2, -18:0.5:-12, AxesLims, 3, 1);
%   imagesc(Anatomy(:,:,1));
%
%   AxesLims = [0, 16; -30, 10; -20, 20];
%   [Anatomy, Structures] = EN_GetAnatomySlices('inia19', 2, -18:0.5:-12, AxesLims, 3, 1);
%   [Anatomy, Structures] = EN_GetAnatomySlices('D99', 2, -18:0.5:-12, AxesLims, 3, 1);
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy © Copyleft 2014-2016, GNU General Public License
%==========================================================================

%================ CHECK INPUTS
if ~all(size(AxesLims)== [3,2])
    error('AxesLims input must be a 3 row by 2 column matrix')
end

switch lower(SubjectID)
    case lower({'D99','Saleem'})
       	Defaults.MRI    = '/Volumes/PROJECTS/murphya/EN_data/Atlases/D99/D99_template.nii'; 
        Defaults.VTKdir = '/Volumes/PROJECTS/murphya/EN_data/Atlases/D99/VTKs';
        Thresh          = 200;
    case lower({'neuromaps','inia19'})
        Defaults.MRI    = '/Volumes/PROJECTS/murphya/EN_data/Atlases/inia19/inia19-t1.nii'; 
        Defaults.VTKdir = '/Volumes/PROJECTS/murphya/EN_data/Atlases/inia19/VTKs';
        Thresh          = 160; 
    case lower('Frey')
      	Defaults.MRI    = '/Volumes/PROJECTS/murphya/EN_data/Atlases/Frey/rhesus_7_model-MNI.nii'; 
        Defaults.VTKdir = '/Volumes/PROJECTS/murphya/EN_data/Atlases/Frey/VTKs';
    case lower('Dexter')
     	Thresh      = 15000;   
        Defaults    = ENT_LoadDefaults(SubjectID);
    case lower('Layla')
        Thresh      = 10000;   
        Defaults    = ENT_LoadDefaults(SubjectID);
    otherwise
        error('Unknown atlas/ subject: %s', SubjectID);
end
StructVolumes   = wildcardsearch(Defaults.VTKdir, '*ROI.nii');


%================ GET REQUESTED SLICES FROM MRI
if ismember(Type, [1,3])
	nii         = load_nii(Defaults.MRI);
    T           = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];
    OriginVox  	= nii.hdr.hist.originator(1:3);
    VoxSize     = nii.hdr.dime.pixdim(2:4);
    VolSize     = nii.hdr.dime.dim(2:4);
 	OriginMM    = OriginVox.*VoxSize;                     
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
             for S = 1:numel(SlicePos)
                B = bwboundaries(Structures{m}(:,:,S));
                 for k = 1:length(B)
                 	switch Plane    %============ Convert coordinates from image pixels to mm from origin
                        case 1
                         	B{k}(:,1) = (B{k}(:,1) + AxesLims(3,1)/VoxSize(3)-1)*VoxSize(3);
                            B{k}(:,2) = (B{k}(:,2) + AxesLims(2,1)/VoxSize(2)-1)*VoxSize(2);
                        case 2
                            B{k}(:,1) = (B{k}(:,1) + AxesLims(3,1)/VoxSize(3)-1)*VoxSize(3);
                            B{k}(:,2) = (B{k}(:,2) + AxesLims(1,1)/VoxSize(1)-1)*VoxSize(1);
                        case 3
                        	B{k}(:,1) = (B{k}(:,1) + AxesLims(2,1)/VoxSize(2)-1)*VoxSize(2);
                            B{k}(:,2) = (B{k}(:,2) + AxesLims(1,1)/VoxSize(1)-1)*VoxSize(1);
                    end
                     
                 	Outlines{S,m,k} = B{k}(:,[2,1]);
                 end
             end
        end
        
    end
    
   
end


%================ PLOT DATA
if TestPlot == 1
    figure('name', sprintf('EN_GetAnatomySlices: %s', SubjectID));
    axh             = tight_subplot(2, ceil(numel(SlicePos)/2), 0.02, 0.02, 0.02);
    StructColors    = jet(numel(Structures));
    MaskAlpha       = 0.3; 
    FillStruct      = 1;
    
    for S = 1:numel(SlicePos)
        axes(axh(S));
        AnatomyRGB{S}   = repmat(double(Anatomy(:,:,S))/double(max(max(Anatomy(:,:,S)))), [1,1,3]);
        imh(S)          = image(AnatomyRGB{S});
        hold on;
        
        %================== Plot structure overlay
        if exist('Structures','var')
            for N = 1:numel(Structures)
                if FillStruct == 1
                    StructBinary{N}(:,:,S) = double(Structures{N}(:,:,S))/double(max(max(Structures{N}(:,:,S))));        	% Normalize mask values
                    imsh(S,N) = imagesc(StructBinary{N}(:,:,S)*N,'alphadata', StructBinary{N}(:,:,S)*MaskAlpha);              % Plot structure filled
                end
                for k = 1:size(Outlines, 3)
                    if ~isempty(Outlines{S,N,k})
                        StructLineH(S,N) = plot(Outlines{S,N,k}(:,1), Outlines{S,N,k}(:,2),'-w','linewidth',1);     % Plot structure boundary outline
                        set(StructLineH(S,N), 'color', StructColors(N,:));
                    end
                end
            end
        end
        
        %================== Scale images to mm
        switch Plane
            case 1
                set(imh(S),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
                set(imsh(S,:),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
            case 2
                set(imh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
                set(imsh(S,:),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
            case 3
                set(imh(S),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
                set(imsh(S,:),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
        end

        axis equal tight xy
        title(sprintf('%.1f mm', SlicePos(S)));
    end
    
    %================== Plot legend
    axes(axh(S+1));
    for S = 1:numel(StructVolumes)
        X = [0,0,1,1];
        Y = [0,1,1,0]-(S*1.5);
        patch(X,Y,ones(size(X))*S);
        [~,Filename] = fileparts(StructVolumes{S});
        Indx = strfind(Filename, '_ROI')-1;
        Filename(strfind(Filename, '_')) = ' ';
        Label{S} = Filename(1:Indx);
        text(X(3)+1, mean(Y([1,2])), Label{S},'fontsize',14);
    end
    set(axh(S+1),'xlim',[-10 10],'ylim',[-10 10]);
    colormap jet;
    axis tight equal off;
end