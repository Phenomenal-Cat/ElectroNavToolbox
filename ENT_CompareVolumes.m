
%========================== ENT_CompareVolumes.m ========================== 
% This simple GUI loads two 3D volumes and allows the user to interactively
% explore and compare them. This can be useful for checking registration of
% volumes or changes between sessions (e.g. pre- vs post-lesion).
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, � Copyleft 2015, GNU General Public License
%==========================================================================

function ENT_CompareVolumes(NiiFile1, NiiFile2)
global Fig Nii

% if nargin < 2
%     [File1, Path] = uigetfile({'*.nii;*.img'},'Select first MR volume');
%     NiiFile1 = fullfile(Path, File1);
%     [File2, Path] = uigetfile({'*.nii;*.img'},'Select second MR volume', Path);
%     NiiFile2 = fullfile(Path, File2);
% else
%   
% end



NiiFile1 = '/Volumes/rawdata/murphya/MRI/Layla/20150515_preElgiloy/r20150519_PreElgiloy_FLASH_025mm_Aligned_BET.nii';
NiiFile2 = '/Volumes/RAWDATA/murphya/MRI/Layla/20150602_post_elgiloy/r20150602_FLASH_postelgiloy_iso_BET_aligned_rot.nii';


%======================= LOAD VOLUMES
Nii{1} = load_nii(NiiFile1);
Nii{2} = load_nii(NiiFile2);

Fig.VoxelDim = Nii{1}.hdr.dime.pixdim(2:4);
Fig.Origin = Nii{1}.hdr.hist.originator(1:3);
Fig.OriginMM = Nii{1}.hdr.hist.originator(1:3).*Fig.VoxelDim;
Fig.VolumeDim = size(Nii{1}.img);
Fig.VolumeDimMM = size(Nii{1}.img).*Fig.VoxelDim;
Fig.AxLims = [-Fig.OriginMM; Fig.VolumeDimMM-Fig.OriginMM]';
Fig.SliceOrientation = 1;
Fig.SlicePosition = Fig.Origin;
Fig.SlicePositionMM = [0 0 0];

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
Fig.MenuLabels = {'Slice orientation','Slice position','Zoom','Pan','Vol1 threshold', 'Vol2 threshold'};
Fig.MenuInputTypes = {'popup','slider','toggle','toggle','slider','slider'};
Fig.MenuInputStrings = {{'Sagittal','Coronal','Axial'}, [], 'on','on',[],[]};
Ypos = 20:25:170;
for n = 1:numel(Fig.MenuLabels)
    Fig.MenuLabelH(n) = uicontrol('Style','text','string',Fig.MenuLabels{n},'pos',[10, Ypos(n), 100, 25],'parent',Fig.UIpannel);
    Fig.MenuInputH(n) = uicontrol('Style',Fig.MenuInputTypes{n},'string',Fig.MenuInputStrings{n}, 'HorizontalAlignment','Left','pos',[120, Ypos(n), 100, 25],'parent',Fig.UIpannel,'Callback',{@MenuInput,n});
end
set(Fig.MenuLabelH, 'HorizontalAlignment','Left', 'backgroundcolor', Fig.Background);
set(Fig.MenuLabelH(2),'min', 1,'max', Fig.VolumeDim(Fig.SliceOrientation),'SliderStep', repmat(1/Fig.VolumeDim(Fig.SliceOrientation),[1,2]), 'value',Fig.SlicePosition(Fig.SliceOrientation));
set(Fig.MenuLabelH(1), 'value', Fig.SliceOrientation);      
set(Fig.MenuInputH([5,6]),'value',1);


%======================= PLOT DATA               
Fig.ax(1) = axes('units','normalized','position',[0.05 0.1 0.3 0.6]);
Fig.Im(1) = imagesc(Fig.AxLims(2,1):Fig.AxLims(2,2), Fig.AxLims(1,1):Fig.AxLims(1,2), Nii{1}.img(:,:,Nii{1}.hdr.hist.originator(3)));
axis equal tight;
hold on;
Fig.OriginHx(1) = plot(xlim, [0 0],'--c');
Fig.OriginHy(1) = plot([0 0], ylim,'--c');
[path, name, ext] = fileparts(NiiFile1);
name(strfind(name,'_'))=' ';
title(name, 'fontsize',18);

Fig.ax(2) = axes('units','normalized','position',[0.4 0.1 0.3 0.6]);
Fig.Im(2) = imagesc(Fig.AxLims(2,1):Fig.AxLims(2,2), Fig.AxLims(1,1):Fig.AxLims(1,2), Nii{2}.img(:,:,Nii{2}.hdr.hist.originator(3)));
axis equal tight;
hold on;
Fig.OriginHx(2) = plot(xlim, [0 0],'--c');
Fig.OriginHy(2) = plot([0 0], ylim,'--c');
[path, name, ext] = fileparts(NiiFile2);
name(strfind(name,'_'))=' ';
title(name, 'fontsize',18);

linkaxes(Fig.ax([1,2]));
colormap gray;
set(Fig.Im, 'ButtonDownFcn',@PointerCallback);

Fig.ax(3) = axes('parent',Fig.UIpannel, 'units','normalized','position',[0.2, 0.6, 0.7, 0.15]);
hist(double(Nii{1}.img(:)), 1000)
set(gca,'ylim',[0 0.005*10^7])

Fig.ax(4) = axes('parent',Fig.UIpannel, 'units','normalized','position',[0.2, 0.4, 0.7, 0.15]);
hist(double(Nii{2}.img(:)), 1000)
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
            set(Fig.MenuLabelH(2),'max',Fig.VolumeDim(Fig.SliceOrientation));
            switch Fig.SliceOrientation
                case 1
                    set(Fig.ax([1,2]), 'xlim', Fig.AxLims(2,:), 'ylim', Fig.AxLims(3,:));
                case 2
                    set(Fig.ax([1,2]), 'xlim', Fig.AxLims(1,:), 'ylim', Fig.AxLims(3,:));
                case 3
                	set(Fig.ax([1,2]), 'xlim', Fig.AxLims(2,:), 'ylim', Fig.AxLims(1,:));
            end
            
        case 2      %========== Change slice position
            Fig.SlicePosition(Fig.SliceOrientation) = round(get(hObj,'Value')*Fig.VolumeDim(Fig.SliceOrientation));
            Fig.SlicePositionMM(Fig.SliceOrientation) = (Fig.SlicePosition(Fig.SliceOrientation).*Fig.VoxelDim(Fig.SliceOrientation))-Fig.Origin(Fig.SliceOrientation);
            set(Fig.MenuLabelH(2), 'string', sprintf('Position = %.2f mm', Fig.SlicePositionMM(Fig.SliceOrientation)));
            UpdateSlice;

        case 3      %========== Zoom
            zoom on;
            
        case 4      %========== Pan
            
        case 5      %========== Adjust intensity thershold
            Thresh(1) = get(hObj,'Value')*max(double(Nii{1}.img(:))); 
            set(Fig.ax(3),'xlim',[0, Thresh(1)]);
            set(Fig.ax(1),'clim',[0, Thresh(1)]);
            
        case 6
          	Thresh(2) = get(hObj,'Value')*max(double(Nii{2}.img(:))); 
            set(Fig.ax(4),'xlim',[0, Thresh(2)]);
            set(Fig.ax(2),'clim',[0, Thresh(2)]);
    end

end

function UpdateSlice
global Fig Nii
    switch Fig.SliceOrientation
        case 1      %======= SAGITTAL
            set(Fig.Im, 'xdata', Fig.AxLims(2,1):Fig.AxLims(2,2), 'ydata', Fig.AxLims(3,1):Fig.AxLims(3,2));
            set(Fig.Im(1), 'cdata', rot90(squeeze(Nii{1}.img(Fig.SlicePosition(Fig.SliceOrientation),:,:))));
            set(Fig.Im(2), 'cdata', rot90(squeeze(Nii{2}.img(Fig.SlicePosition(Fig.SliceOrientation),:,:))));
            
            
        case 2      %======= CORONAL
            set(Fig.Im, 'xdata', Fig.AxLims(1,1):Fig.AxLims(1,2), 'ydata', Fig.AxLims(3,1):Fig.AxLims(3,2));
            set(Fig.Im(1), 'cdata', rot90(squeeze(Nii{1}.img(:,Fig.SlicePosition(Fig.SliceOrientation),:))));
            set(Fig.Im(2), 'cdata', rot90(squeeze(Nii{2}.img(:,Fig.SlicePosition(Fig.SliceOrientation),:))));
            
            
        case 3      %======= AXIAL
          	set(Fig.Im, 'xdata', Fig.AxLims(2,1):Fig.AxLims(2,2), 'ydata', Fig.AxLims(1,1):Fig.AxLims(1,2));
            set(Fig.Im(1), 'cdata', Nii{1}.img(:,:,Fig.SlicePosition(Fig.SliceOrientation)));
            set(Fig.Im(2), 'cdata', Nii{2}.img(:,:,Fig.SlicePosition(Fig.SliceOrientation)));
    end
    
    %======= CURRENT VOXEL DATA
    
    
    
end

function PointerCallback(objectHandle, eventData)
global Fig
	axesHandle  = get(objectHandle,'Parent');
    coordinates = get(axesHandle,'CurrentPoint');
    CoordinateMM = coordinates(1,:);
%     CoordinateVox = CoordinateMM-Fig.*Fig.VoxelDim;
    
    if ~isfield(Fig, 'SelectedVoxH')
        for a = 1:2
            axes(Fig.ax(a));
            Fig.SelectedVoxH(1,a) = plot(repmat(CoordinateMM(1), [1,2]), Fig.AxLims(3,:), '-r');
            Fig.SelectedVoxH(2,a) = plot(Fig.AxLims(1,:), repmat(CoordinateMM(2), [1,2]), '-r');
        end
        set(Fig.SelectedVoxH, 'PickableParts','none');
    else
        set(Fig.SelectedVoxH(1,:),'xdata',repmat(CoordinateMM(1), [1,2]));
        set(Fig.SelectedVoxH(2,:),'ydata',repmat(CoordinateMM(2), [1,2]));
    end
    
end
