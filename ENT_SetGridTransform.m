
%====================== ENT_SetGridTransform.m ============================
% This function loads and displays an anatomical MRI with gadolinium-filled
% grid in place. The user may then manually adjust the position and
% orientation of the grid to match the scan, and the resulting
% transformation matrix is saved.
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function ENT_SetGridTransform(NiiFile)


clear all
global Fig Grid nii FV

%============= LOAD DATA
if nargin == 0
%     [file, path]        = uigetfile({'*.nii'}, 'Select ACPC-aligned grid scan');
%     NiiFile             = fullfile(path, file);
%     NiiFile             = '/Volumes/rawdata/murphya/MRI/Dexter/20151203/Dexter_20151203_MDEFT_hi_ACPC.nii';
    NiiFile             = '/Volumes/rawdata/murphya/MRI/Dexter/20151209/Dexter_20151209_ACPC.nii';
end
nii = LoadMRI(NiiFile);


Grid.VolFile           	= 'Grid1.nii';
Grid.SurfFile         	= 'Grid1.stl';
Grid.nii              	= load_nii(Grid.VolFile);                 	% Load grid volume
[FV.vertices, FV.faces]	= stlread(Grid.SurfFile);                  	% Load grid surface
Grid.rot                = [-15, 0, 0];                              % Rotations about cardinal axes (degrees)
Grid.Trans              = [7, -15, 15];                             % Translations relative to AC origin
Grid.Tform             	= makehgtform('translate',Grid.Trans(1), Grid.Trans(2), Grid.Trans(3),'xrotate',Grid.rot(1),'yrotate',Grid.rot(2),'zrotate',Grid.rot(3));
FV(2).vertices          = ENT_ApplyTform(Grid.Tform, FV(1).vertices);
Grid.Colors            	= [1 0 0; 0 1 0];


%% ========================== OPEN FIGURE WINDOW ==========================
Fig.Handle      = figure('position', get(0, 'ScreenSize'),'name', 'ENT_SetGridTransform');
Fig.axh         = tight_subplot(2,3,0.05, 0.05, 0.05);
delete(Fig.axh([3,6]));
Fig.axh([3,6])  = [];
Fig.UIFontsize  = 14;
Fig.Background  = get(Fig.Handle, 'color');


%=============== MRI GUI PANEL
Fig.MRI.BoxPos      = [0.75,0.6,0.2,0.3];
Fig.MRI.InputDim    = [150 20];
Fig.MRI.Labels      = {'Anatomical','Resolution (mm)','Origin (vx)','Permute','Flip', 'Rotate 90', 'Contrast'};
Fig.MRI.Style       = {'pushbutton','text','edit','pushbutton','pushbutton','pushbutton','pushbutton'};
Fig.MRI.List        = {nii.filename, nii.VoxDim, nii.OriginVox, 'Permute', 'Flip', 'Rotate 90', 'Contrast'};
Fig.MRI.UIhandle    = uipanel('Title','Anatomical MRI','FontSize', Fig.UIFontsize,'BackgroundColor',Fig.Background,'Units','normalized','Position',Fig.MRI.BoxPos);
for i = 1:numel(Fig.MRI.Labels)
    Pos = numel(Fig.MRI.Labels)-i;
    Fig.MRI.LabelPos{i} = [20, 10+Pos*Fig.MRI.InputDim(2)*1.5, Fig.MRI.InputDim];
    Fig.MRI.LabelHandle(i) = uicontrol( 'Style','Text',...
                                        'String',Fig.MRI.Labels{i},...
                                        'HorizontalAlignment','Left',...
                                        'pos',Fig.MRI.LabelPos{i},...
                                        'parent',Fig.MRI.UIhandle);
   	if ~any(ismember(i, [2,3]))
        Fig.MRI.InputHandle(i) = uicontrol( 'Style',Fig.MRI.Style{i},...
                                            'String',Fig.MRI.List{i},...
                                            'HorizontalAlignment','Left',...
                                            'pos',[Fig.MRI.InputDim(1)+10,15+Pos*Fig.MRI.InputDim(2)*1.5,100,20],...
                                            'Callback',{@MRIparams,i,0},...
                                            'parent',Fig.MRI.UIhandle);
    elseif any(ismember(i, [2,3]))
        for n = 1:3
            Fig.MRI.MultiInputH(i,n) = uicontrol( 'Style',Fig.MRI.Style{i},...
                                            'String', num2str(Fig.MRI.List{i}(n)),...
                                            'HorizontalAlignment','Left',...
                                            'pos',[Fig.MRI.InputDim(1)+10+((n-1)*50), 15+Pos*Fig.MRI.InputDim(2)*1.5, 40, 20],...
                                            'HorizontalAlignment','Left',...
                                            'Callback',{@MRIparams,i,n},...
                                            'parent',Fig.MRI.UIhandle);
        end
    end
end
% set([Fig.MRI.LabelHandle, Fig.MRI.InputHandle], 'BackgroundColor',Fig.Background);
% set(Fig.MRI.InputHandle(2),'value', Fig.MRI.CurrentDate);
% set(Fig.MRI.InputHandle(3),'value', Electrode(1).Selected);
% set(Fig.MRI.InputHandle(4),'value', find(~cellfun(@isempty, strfind(Electrode(Electrode(1).Selected).AllTypes, Electrode(Electrode(1).Selected).Brand))));

%=============== GRID GUI PANEL
Fig.Grid.BoxPos      = [0.75,0.2,0.2,0.3];
Fig.Grid.InputDim    = [150 20];
Fig.Grid.Labels      = {'Grid type','Origin (mm)','Rotation','Transparency','Save transform'};
Fig.Grid.Style       = {'popupmenu','edit','edit','edit','pushbutton'};
Fig.Grid.List        = {{'19mm cylindrical'}, Grid.Trans, Grid.rot, 1, 'Save'};
Fig.Grid.UIhandle    = uipanel('Title','Grid','FontSize', Fig.UIFontsize,'BackgroundColor',Fig.Background,'Units','normalized','Position',Fig.Grid.BoxPos);
for i = 1:numel(Fig.Grid.Labels)
    Pos = numel(Fig.Grid.Labels)-i;
    Fig.Grid.LabelPos{i} = [20, 10+Pos*Fig.Grid.InputDim(2)*1.5, Fig.Grid.InputDim];
    Fig.Grid.LabelHandle(i) = uicontrol( 'Style','Text',...
                                        'String',Fig.Grid.Labels{i},...
                                        'HorizontalAlignment','Left',...
                                        'pos',Fig.Grid.LabelPos{i},...
                                        'parent',Fig.Grid.UIhandle);
    if ~any(ismember(i, [2,3]))                                
        Fig.Grid.InputHandle(i) = uicontrol( 'Style',Fig.Grid.Style{i},...
                                            'String',Fig.Grid.List{i},...
                                            'HorizontalAlignment','Left',...
                                            'pos',[Fig.Grid.InputDim(1)+10,15+Pos*Fig.Grid.InputDim(2)*1.5,100,20],...
                                            'Callback',{@Gridparams,i},...
                                            'parent',Fig.Grid.UIhandle);
    elseif any(ismember(i, [2,3]))
        for n = 1:3
            Fig.Grid.MultiInputH(i,n) = uicontrol( 'Style',Fig.Grid.Style{i},...
                                        'String',Fig.Grid.List{i}(n),...
                                        'HorizontalAlignment','Left',...
                                        'pos',[Fig.Grid.InputDim(1)+10+((n-1)*50),15+Pos*Fig.Grid.InputDim(2)*1.5,40,20],...
                                        'Callback',{@Gridparams,i},...
                                        'parent',Fig.Grid.UIhandle);
        end
    end
end
set([Fig.Grid.LabelHandle,Fig.Grid.InputHandle([1,4])], 'BackgroundColor',Fig.Background);



UpdatePlots;


end

%% ========================== SUBFUNCTIONS ================================

%============= Load MRI volume
function nii = LoadMRI(NiiFile)
    nii                 = load_nii(NiiFile);                            % Load anatomical volume
    [~, nii.filename]   = fileparts(nii.fileprefix);
    nii.OriginVox     	= nii.hdr.hist.originator(1:3); 
    nii.VoxDim          = nii.hdr.dime.pixdim(2:4);
    nii.DimVox          = nii.hdr.dime.dim(2:4);
    nii.DimMm           = nii.VoxDim.*nii.DimVox;
    nii.OriginMm        = nii.OriginVox.*nii.VoxDim;
    nii.AxLims          = [-nii.OriginMm; nii.DimMm-nii.OriginMm];
    nii.sform           = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z];
end

%============= Plot slices
function UpdatePlots

global Fig Grid nii FV

%     ph(1) = patch(nii.AxLims([1,1,2,2],1)', nii.AxLims([1,2,2,1],2)', zeros(1,4), zeros(1,4));%, 'Cdata', squeeze(nii.img(:,:,nii.OriginVox(3))),'FaceColor','texturemap');
%     ph(2) = patch(zeros(1,4), nii.AxLims([1,1,2,2],2)', nii.AxLims([1,2,2,1],3)', zeros(1,4));%, 'Cdata', squeeze(nii.img(nii.OriginVox(1),:,:)),'FaceColor','texturemap');
%     ph(3) = patch(nii.AxLims([1,1,2,2],1)', zeros(1,4), nii.AxLims([1,2,2,1],3)', zeros(1,4));%, 'Cdata', squeeze(nii.img(:,nii.OriginVox(2),:)),'FaceColor','texturemap');

    %============= PLOT 3D VIEW
    axes(Fig.axh(1));
    for v = 1:3
        switch v                                             
            case 1
                XCoords = repmat(nii.AxLims(1,1), [2,2]);
                YCoords = [nii.AxLims(1,2), nii.AxLims(1,2); nii.AxLims(2,2), nii.AxLims(2,2)]';
                ZCoords = [nii.AxLims(:,3), nii.AxLims(:,3)];
                CurrentSlice = repmat(squeeze(nii.img(:,:,nii.OriginVox(3))),[1,1,3]);
            case 2
                XCoords = [nii.AxLims(1,1), nii.AxLims(1,1); nii.AxLims(2,1), nii.AxLims(2,1)]';
                YCoords = repmat(nii.AxLims(1,2), [2,2]);
                ZCoords = [nii.AxLims(:,3), nii.AxLims(:,3)];
                CurrentSlice = repmat(squeeze(nii.img(nii.OriginVox(1),:,:)),[1,1,3]);
            case 3

                XCoords = [nii.AxLims(1,1), nii.AxLims(1,1); nii.AxLims(2,1), nii.AxLims(2,1)]';
                YCoords = [nii.AxLims(:,2), nii.AxLims(:,2)];
                ZCoords = repmat(nii.AxLims(1,3), [2,2]);
                CurrentSlice = repmat(squeeze(nii.img(:,nii.OriginVox(2),:)),[1,1,3]);
        end    
        Fig.ph(v) = surf(XCoords, YCoords, ZCoords,'CData', repmat(squeeze(nii.img(:,:,nii.OriginVox(3))),[1,1,3]),'FaceColor','texturemap','EdgeColor','k');         % Draw MRI slice to axes
    end
    
    %=========== PLOT 3D GRID(S)
    for g = 1:2
        Grid.h(g) = patch(FV(g), 'facecolor', Grid.Colors(g,:), 'edgecolor', 'none');     % Plot grid at origin
        hold on;
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
    Fig.lineh(2,1) = plot([0 0], ylim, '--w');
    Fig.lineh(2,2) = plot(xlim, [0 0], '--w');
    xlabel('Y (mm)', 'fontsize', 14);
    ylabel('X (mm)', 'fontsize', 14);
    title('Axial', 'fontsize', 16);
    axis equal tight xy;

    axes(Fig.axh(3));
    Fig.imh(3) = imagesc(nii.AxLims(:,2)', nii.AxLims(:,3)', squeeze(nii.img(nii.OriginVox(1),:,:))');
    hold on;
    Fig.lineh(3,1) = plot([0 0], ylim, '--w');
    Fig.lineh(3,2) = plot(xlim, [0 0], '--w');
    xlabel('Y (mm)', 'fontsize', 14);
    ylabel('Z (mm)', 'fontsize', 14);
    title('Sagittal', 'fontsize', 16);
    axis equal tight xy;

    axes(Fig.axh(4));
    Fig.imh(4) = imagesc(nii.AxLims(:,1)', nii.AxLims(:,3)', squeeze(nii.img(:,nii.OriginVox(2),:))');
    hold on;
    Fig.lineh(4,1) = plot([0 0], ylim, '--w');
    Fig.lineh(4,2) = plot(xlim, [0 0], '--w');
    xlabel('X (mm)', 'fontsize', 14);
    ylabel('Z (mm)', 'fontsize', 14);
    title('Coronal', 'fontsize', 16);
    axis equal tight xy;

    colormap gray
end


function UpdateSlices(XYZ)
global Fig Grid nii
    set(Fig.imh(2), 'cdata', squeeze(nii.img(:,:,XYZ(3))), 'xdata', nii.AxLims(:,2)', 'ydata', nii.AxLims(:,1)');
    set(Fig.imh(3), 'cdata', squeeze(nii.img(XYZ(1),:,:))', 'xdata', nii.AxLims(:,2)', 'ydata', nii.AxLims(:,3)');
    set(Fig.imh(4), 'cdata', squeeze(nii.img(:,XYZ(2),:))', 'xdata', nii.AxLims(:,1)', 'ydata', nii.AxLims(:,3)');

    set(Fig.lineh(2,1), 'ydata', nii.AxLims([1,2],1)');
    set(Fig.lineh(2,2), 'xdata', nii.AxLims([1,2],2)');
 	set(Fig.lineh(3,1), 'ydata', nii.AxLims([1,2],3)');
    set(Fig.lineh(3,2), 'xdata', nii.AxLims([1,2],2)');
 	set(Fig.lineh(4,1), 'ydata', nii.AxLims([1,2],3)');
    set(Fig.lineh(4,2), 'xdata', nii.AxLims([1,2],1)');
    
    for ax = 2:4
        axes(Fig.axh(ax));
        axis equal tight
    end
end


function MRIparams(hObj, event, indx, indx2)
global Fig Grid nii

switch indx
    case 1          %=============== LOAD NEW MRI VOLUME
        [file, path]    = uigetfile('*.nii', 'Select MRI to load');
        NiiFile         = fullfile(path, file);
        nii             = LoadMRI(NiiFile);
        nii.Filename    = file;
        set(Fig.MRI.InputHandle(1), 'string', nii.Filename);
     	for n = 1:3
        	set(Fig.MRI.MultiInputH(2,n),'String', num2str(nii.VoxDim(n)));
            set(Fig.MRI.MultiInputH(3,n),'String', num2str(nii.OriginVox(n)));
        end
        UpdateSlices(nii.OriginVox);
        
    case 3          %=============== ADJUST ORIGIN COORDINATES
        nii.OriginVox(indx2)= str2num(get(hObj,'string'));
        nii.OriginMm        = nii.OriginVox.*nii.VoxDim;
        nii.AxLims          = [-nii.OriginMm; nii.DimMm-nii.OriginMm];
        UpdateSlices(nii.OriginVox);

        
    case 4          %===============           
        
        
end


end

function Gridparams(hObj, event, indx)
global Fig Grid nii

switch indx
    case 1          %=============== CHANGE GRID TYPE
        
    case {2,3,4}    %=============== CHANGE ORIGIN
        
    case {5,6,7}    %=============== CHANGE ORIENTATION
        
    case 8          %=============== CHANGE OPACITY
        
end
end