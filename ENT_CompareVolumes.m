
%========================== ENT_CompareVolumes.m ========================== 
% This simple GUI loads two 3D volumes and allows the user to interactively
% explore and compare them. This can be useful for checking registration of
% volumes or changes between sessions (e.g. pre- vs post-lesion).
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function ENT_CompareVolumes(NiiFile1, NiiFile2)
global Fig Nii

if nargin < 2
    DefaultPath = '/Volumes/rawdata/murphya/MRI/';
    [File1, Path] = uigetfile({'*.nii;*.img'},'Select first MR volume', DefaultPath);
    NiiFile1 = fullfile(Path, File1);
    [File2, Path] = uigetfile({'*.nii;*.img'},'Select second MR volume', Path);
    NiiFile2 = fullfile(Path, File2);
end
NiiFile{1} = NiiFile1;
NiiFile{2} = NiiFile2;

% % NiiFile1 = '/Volumes/rawdata/murphya/MRI/Layla/20150515_preElgiloy/r20150519_PreElgiloy_MDEFT_025mm_BET.nii';
% % NiiFile2 = '/Volumes/RAWDATA/murphya/MRI/Layla/20150602_post_elgiloy/r20150602_MDEFT_postelgiloy_025mm_BET_Masked_rot.nii';
% 
% % NiiFile1 = '/Volumes/rawdata/murphya/MRI/Layla/20150515_preElgiloy/r20150519_PreElgiloy_FLASH_025mm_Aligned_BET.nii';
% % NiiFile2 = '/Volumes/RAWDATA/murphya/MRI/Layla/20150602_post_elgiloy/r20150602_FLASH_postelgiloy_iso_BET_aligned_rot.nii';
% 
% % NiiFile1 = '/Volumes/rawdata/murphya/MRI/Dexter/20150917/Dexter_20150917_ACPC_resliced.nii';
% NiiFile2 = '/Volumes/rawdata/murphya/MRI/Dexter/20151209/Dexter_20151209_ACPC.nii';
% NiiFile1 = '/Volumes/rawdata/murphya/MRI/Dexter/20151209/Dexter_20151209_ACPC.nii';
% % NiiFile2 = '/Volumes/RAWDATA/murphya/MRI/Dexter/20140219/Dexter_meanT1_BET_resliced.nii';


%======================= LOAD VOLUMES
for n = 1:2
    fprintf('Loading volume %d: %s...\n', n, NiiFile{n})
    Nii{n} = load_nii(NiiFile{n});
end

if Nii{1}.hdr.dime.pixdim(2:4) ~= Nii{2}.hdr.dime.pixdim(2:4)
    fprintf('WARNING: Volume resolution mismatch!\n')
    for n = 1:2
       fprintf('\t%s voxel size = [%.2f, %.2f, %.2f] mm\n', NiiFile{n}, Nii{n}.hdr.dime.pixdim(2:4));
    end
    [~,FileToReslice] = max([prod(Nii{1}.hdr.dime.pixdim(2:4)), prod(Nii{2}.hdr.dime.pixdim(2:4))]);
    FileNotToReslice = find([1,2]~= FileToReslice);
	fprintf('Attempting to reslice lower res volume (%s) to match higher res...\n', NiiFile{FileToReslice});
    [path,file,ext] = fileparts(NiiFile{FileToReslice});
    ReslicedFilename = fullfile(path, sprintf('%s_resliced%s', file, ext));
    Verbose = 1;        % Show progress of reslice
    Method  = 1;        % Trilinear interpolation
	reslice_nii(NiiFile{FileToReslice}, ReslicedFilename, Nii{FileNotToReslice}.hdr.dime.pixdim(2:4), Verbose, Method);
    NiiFile{FileToReslice} = ReslicedFilename;
    Nii{FileToReslice} = load_nii(NiiFile{FileToReslice});
end

for n = 1:2
    Fig.VoxelDim{n}       	= Nii{n}.hdr.dime.pixdim(2:4);
    Fig.Origin{n}          	= round(Nii{n}.hdr.hist.originator(1:3));
    Fig.OriginMM{n}        	= Nii{n}.hdr.hist.originator(1:3).*Fig.VoxelDim{n};
    Fig.VolumeDim{n}       	= size(Nii{n}.img);
    Fig.VolumeDimMM{n}    	= size(Nii{n}.img).*Fig.VoxelDim{n};
    Fig.AxLims{n}          	= [-Fig.OriginMM{n}; Fig.VolumeDimMM{n}-Fig.OriginMM{n}]';
    Fig.SliceOrientation    = 1;
    Fig.SlicePosition       = Fig.Origin{1};
    Fig.SlicePositionMM    	= [0 0 0];
end

%======================= OPEN FIGURE
Fig.Background = [0.5,0.5,0.5];
Fig.scnsize = get(0,'ScreenSize');                                              % Get screen resolution
Fig.Rect = [0 0 Fig.scnsize(3), Fig.scnsize(4)];                                % Set figure winow to fullscreen
Fig.Handle = figure('Name','Compare Volumes',...                                % Open a figure window with specified title
                    'Color',Fig.Background,...                                  % Set the figure window background color
                    'Renderer','OpenGL',...                                     % Use OpenGL renderer
                    'Position', Fig.Rect,...                                    % position figure window to fit fullscreen
                    'NumberTitle','off',...                                     % Remove figure number from title
                    'IntegerHandle','off');                                     % Don't use integer handles
Fig.UIpannel = uipanel('Title','Menu','FontSize',14,'BackgroundColor',Fig.Background,'Units','normalized','Position',[0.75, 0.2, 0.2, 0.6],'Parent',Fig.Handle);
Fig.MenuLabels = {'Slice orientation','Slice position','Rotate 90','Flip','Vol1 threshold', 'Vol2 threshold'};
Fig.MenuInputTypes = {'popup','slider','pushbutton','pushbutton','slider','slider'};
Fig.MenuInputStrings = {{'Sagittal','Coronal','Axial'}, [], 'on','on',[],[]};
Ypos = 20:25:195;
for n = 1:numel(Fig.MenuLabels)
    Fig.MenuLabelH(n) = uicontrol('Style','text','string',Fig.MenuLabels{n},'pos',[10, Ypos(n), 100, 25],'parent',Fig.UIpannel);
    Fig.MenuInputH(n) = uicontrol('Style',Fig.MenuInputTypes{n},'string',Fig.MenuInputStrings{n}, 'HorizontalAlignment','Left','pos',[120, Ypos(n), 100, 25],'parent',Fig.UIpannel,'Callback',{@MenuInput,n});
end
for n = (1:2)+numel(Fig.MenuLabels)
    Fig.MenuInputH(n) = uicontrol('Style','slider','HorizontalAlignment','Left','pos',[240, Ypos(n-2), 100, 25],'parent',Fig.UIpannel,'Callback',{@MenuInput,n});
end
set(Fig.MenuLabelH, 'HorizontalAlignment','Left', 'backgroundcolor', Fig.Background);
set(Fig.MenuInputH(2),'min', 1,'max', Fig.VolumeDim{1}(Fig.SliceOrientation),'SliderStep', repmat(1/Fig.VolumeDim{1}(Fig.SliceOrientation),[1,2]), 'value',Fig.SlicePosition(Fig.SliceOrientation));
set(Fig.MenuLabelH(1), 'value', Fig.SliceOrientation);      
set(Fig.MenuInputH([5,6]),'value',0);
set(Fig.MenuInputH([7,8]),'value',1);


%======================= PLOT DATA          
Fig.AxPos = {[0.05 0.1 0.3 0.6], [0.4 0.1 0.3 0.6]};
for n = 1:2
    Fig.ax(n) = axes('units','normalized','position',Fig.AxPos{n});
    SliceImage = squeeze(Nii{n}.img(Fig.Origin{n}(Fig.SliceOrientation),:,:));
    Fig.Im(n) = imagesc(Fig.AxLims{n}(2,1):Fig.AxLims{n}(2,2), Fig.AxLims{n}(1,1):Fig.AxLims{n}(1,2), SliceImage);
    axis equal tight xy;
    hold on;
    Fig.OriginHx(n) = plot(xlim, [0 0],'--c');
    Fig.OriginHy(n) = plot([0 0], ylim,'--c');
    [path, name, ext] = fileparts(NiiFile{n});
    name(strfind(name,'_'))=' ';
    title(name, 'fontsize',18);
end

% Fig.ax(2) = axes('units','normalized','position',[0.4 0.1 0.3 0.6]);
% Fig.Im(2) = imagesc(Fig.AxLims{2}(2,1):Fig.AxLims{2}(2,2), Fig.AxLims{2}(1,1):Fig.AxLims{2}(1,2), Nii{2}.img(:,:,Fig.Origin{2}(Fig.SliceOrientation)));
% axis equal tight xy;
% hold on;
% Fig.OriginHx(2) = plot(xlim, [0 0],'--c');
% Fig.OriginHy(2) = plot([0 0], ylim,'--c');
% [path, name, ext] = fileparts(NiiFile{2});
% name(strfind(name,'_'))=' ';
% title(name, 'fontsize',18);

linkaxes(Fig.ax([1,2]));
colormap gray;
set(Fig.Im, 'ButtonDownFcn',@PointerCallback);

Fig.ax(3) = axes('parent',Fig.UIpannel, 'units','normalized','position',[0.2, 0.6, 0.7, 0.15]);
hist(double(Nii{1}.img(:)), 1000)
hold on;
set(gca,'ylim',[0 0.005*10^7]);
Fig.Thresh = repmat(get(gca, 'xlim'),[2,1]);

Fig.ax(4) = axes('parent',Fig.UIpannel, 'units','normalized','position',[0.2, 0.4, 0.7, 0.15]);
hist(double(Nii{2}.img(:)), 1000);
hold on;
set(gca,'ylim',[0 0.005*10^7])
xlabel('Voxel intensity');
ylabel('# voxels');


end


function MenuInput(hObj, Event, Indx)
global Fig Nii
    switch Indx
        case 1      %========== Change orientation
            Fig.SliceOrientation = get(hObj,'Value'); 
            UpdateSlice;
            set(Fig.MenuInputH(2),'max', Fig.VolumeDim(Fig.SliceOrientation),'SliderStep', repmat(1/Fig.VolumeDim(Fig.SliceOrientation),[1,2]), 'value',Fig.SlicePosition(Fig.SliceOrientation));
            switch Fig.SliceOrientation
                case 1

                    for i=1:2
                        axes(Fig.ax(i));
                       	set(Fig.ax(i), 'xlim', Fig.AxLims{i}(2,:), 'ylim', Fig.AxLims{i}(3,:));
                        xlabel('Posterior-Anterior (mm)','fontsize',16);
                        ylabel('Inferior-Superior (mm)','fontsize',16);
                    end
                case 2
                 	for i=1:2
                        axes(Fig.ax(i));
                        set(Fig.ax(i), 'xlim', Fig.AxLims{i}(1,:), 'ylim', Fig.AxLims{i}(3,:));
                        xlabel('Medial-Lateral (mm)','fontsize',16);
                        ylabel('Inferior-Superior (mm)','fontsize',16);
                    end
                case 3
                  	for i=1:2
                        axes(Fig.ax(i));
                        set(Fig.ax(i), 'xlim', Fig.AxLims{i}(2,:), 'ylim', Fig.AxLims{i}(1,:));
                        xlabel('Medial-Lateral (mm)','fontsize',16);
                        ylabel('Posterior-Anterior (mm)','fontsize',16);
                    end
            end
            
            
        case 2      %========== Change slice position
            Fig.SlicePosition(Fig.SliceOrientation) = round(get(hObj,'Value'));
            n = 1;
            Fig.SlicePositionMM = (Fig.SlicePosition-Fig.Origin{n}).*Fig.VoxelDim{n};
            set(Fig.MenuLabelH(2), 'string', sprintf('Position = %.2f mm', Fig.SlicePositionMM(Fig.SliceOrientation)));
            UpdateSlice;

        case 3      %========== Rotate image 90 degrees
            disp('rotating...')
            for im = 1:2
                switch Fig.SliceOrientation
                    case 1
                        Nii{im}.img = permute(Nii{im}.img,[2,3,1]);
                        Nii{im}.img = rot90(Nii{im}.img);
                        Nii{im}.img = permute(Nii{im}.img,[3,1,2]); 
                    case 2
                        Nii{im}.img = permute(Nii{im}.img,[1,3,2]);
                        Nii{im}.img = rot90(Nii{im}.img);
                        Nii{im}.img = permute(Nii{im}.img,[1,3,2]);
                    case 3
                        Nii{im}.img = rot90(Nii{im}.img);
                end
            end
            UpdateSlice;
            
        case 4      %========== Flip image
            for im = 1:2
                switch Fig.SliceOrientation
                    case 1
                        Nii{im}.img = Nii{im}.img(end:-1:1,:,:);

                    case 2
                        Nii{im}.img = Nii{im}.img(:,end:-1:1,:);
                        
                    case 3
                        Nii{im}.img = Nii{im}.img(:,:,end:-1:1);
                end
            end
            
            
        case 5      %========== Vol 1 Lower - Adjust intensity thershold
            Fig.Thresh(1,1) = min(double(Nii{1}.img(:)))+(get(hObj,'Value')*max(double(Nii{1}.img(:))));
            set(Fig.ax(3),'xlim',Fig.Thresh(1,:));
            set(Fig.ax(1),'clim',Fig.Thresh(1,:));
            
        case 6      %========== Vol 2 Lower 
          	Fig.Thresh(2,1) = min(double(Nii{2}.img(:)))+(get(hObj,'Value')*max(double(Nii{2}.img(:))));
            set(Fig.ax(4),'xlim',Fig.Thresh(2,:));
            set(Fig.ax(2),'clim',Fig.Thresh(2,:));
            
        case 7      %========== Vol 1 Upper 
            Fig.Thresh(1,2) = get(hObj,'Value')*max(double(Nii{1}.img(:)));
            set(Fig.ax(3),'xlim',Fig.Thresh(1,:));
            set(Fig.ax(1),'clim',Fig.Thresh(1,:));
            
        case 8      %========== Vol 2 Upper 
           	Fig.Thresh(2,2) = get(hObj,'Value')*max(double(Nii{2}.img(:)));
            set(Fig.ax(4),'xlim',Fig.Thresh(2,:));
            set(Fig.ax(2),'clim',Fig.Thresh(2,:));
    end

end

function UpdateSlice
global Fig Nii
    for n = 1:2
        switch Fig.SliceOrientation
            case 1      %======= SAGITTAL
                set(Fig.Im(n), 'xdata', Fig.AxLims{n}(2,1):Fig.AxLims{n}(2,2), 'ydata', Fig.AxLims{n}(3,1):Fig.AxLims{n}(3,2));
                set(Fig.Im(n), 'cdata', rot90(squeeze(Nii{n}.img(Fig.SlicePosition(Fig.SliceOrientation),:,:))));

            case 2      %======= CORONAL
                set(Fig.Im(n), 'xdata', Fig.AxLims{n}(1,1):Fig.AxLims{n}(1,2), 'ydata', Fig.AxLims{n}(3,1):Fig.AxLims{n}(3,2));
                set(Fig.Im(n), 'cdata', rot90(squeeze(Nii{n}.img(:,Fig.SlicePosition(Fig.SliceOrientation),:))));

            case 3      %======= AXIAL
                set(Fig.Im(n), 'xdata', Fig.AxLims{n}(2,1):Fig.AxLims{n}(2,2), 'ydata', Fig.AxLims{n}(1,1):Fig.AxLims{n}(1,2));
                set(Fig.Im(n), 'cdata', Nii{n}.img(:,:,Fig.SlicePosition(Fig.SliceOrientation)));
        end    
    end
end

function PointerCallback(objectHandle, eventData)
global Fig Nii
	axesHandle  = get(objectHandle,'Parent');
    coordinates = get(axesHandle,'CurrentPoint');
    CoordinateMM = [coordinates(1,[1,2]), Fig.SlicePositionMM(Fig.SliceOrientation)];
    switch Fig.SliceOrientation
        case 2
            CoordinateMM = CoordinateMM([1,3,2]);
        case 3
            CoordinateMM = CoordinateMM([1,3,2]);
    end

    
    %============= GET INTENSITY VALUES OF SELECTED VOXEL
    for im = 1:2
        CoordinateVox{im} = round(CoordinateMM./Fig.VoxelDim{im})+Fig.Origin{im};
        Intensity(im) = Nii{im}.img(CoordinateVox{im}(3), CoordinateVox{im}(2), CoordinateVox{im}(1));
        if ~isfield(Fig, 'IntensityH') || numel(Fig.IntensityH)<2 
            axes(Fig.ax(2+im));
            Fig.IntensityH(im) = plot(repmat(Intensity(im),[1,2]), ylim, '-r');
        else
            set(Fig.IntensityH(im), 'xdata', repmat(Intensity(im),[1,2]));
        end
    end

    %============= UPDATE CROSSHAIRS
    if ~isfield(Fig, 'SelectedVoxH')
        for a = 1:2
            axes(Fig.ax(a));
            Fig.SelectedVoxH(1,a) = plot(repmat(CoordinateMM(1), [1,2]), Fig.AxLims{a}(3,:), '-r');
            Fig.SelectedVoxH(2,a) = plot(Fig.AxLims{a}(1,:), repmat(CoordinateMM(2), [1,2]), '-r');
        end
        set(Fig.SelectedVoxH, 'PickableParts','none');
    else
        set(Fig.SelectedVoxH(1,:),'xdata',repmat(coordinates(1,1), [1,2]));
        set(Fig.SelectedVoxH(2,:),'ydata',repmat(coordinates(1,2), [1,2]));
    end
    
	%======= CURRENT VOXEL DATA
%     Vox.Pos = 
%     VoxelIntensity = 
    
    
end
