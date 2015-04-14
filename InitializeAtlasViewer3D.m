%========================= InitializeAtlasViewer3D ========================
% This function provides a graphical user interface for setting parameters 
% related to 3D scene rendering. Parameters can be saved and loaded, and
% the updated parameters are returned in the structure 'Params'.
%
% INPUTS:
%   DefaultInputs:          optional string containing full path of .mat
%                           file containing previously saved parameters.
% OUTPUT:
%   Params.Stereomode:      Method of stereoscopic presentation (1-10). For
%                           monocular rendering select 0.
%   Params.ViewingDist:     Viewing distance (eyes to screen) in centimetres.
%   Params.IPD:             Interpupillary distance in centimetres.
%   Params.ScreenDims:      2 element vector containing physical screen
%                           dimensions [width, height] in centimetres.
%   Params.ClippingPlanes:  2 element vector containing distance of frustum
%                           clipping planes along the z-axis relative to the
%                           plane of the screen [-near, far] in centimetres.
%   Params.Perspective:     0 = Orthographic, 1 = Perspective projection.
%   Params.LightingOn:      Polygon face lighting on/ off
%   Params.Ambient:         4 element [RGBA] vector specifying ambient lighting
%   Params.Diffuse:         4 element [RGBA] vector specifying diffuse lighting
%   Params.Specular:        4 element [RGBA] vector specifying specular lighting
%   Params.Background:      4 element [RGBA] vector specifying background color
%   Params.
%
%
% REVISIONS:
%   05/04/2014 - Written by APM
%     ___  ______  __   __
%    /   ||  __  \|  \ |  \    APM SUBFUNCTIONS
%   / /| || |__/ /|   \|   \   Aidan P. Murphy - murphyap@mail.nih.gov
%  / __  ||  ___/ | |\   |\ \  Section of Cognitive Neurophysiology and Imaging
% /_/  |_||_|     |_| \__| \_\ Laboratory of Neuropsychology, NIMH
%==========================================================================

function ParamsOut = InitializeAtlasViewer3D(DefaultInputs)

persistent Params Fig;
addpath(genpath('ENSubfunctions'));


if nargin == 0
    Params.Stereomode       = 6;
    Params.ViewingDist      = 50;
    Params.IPD              = 6.4;
    Params.ScreenDims       = [22, 32];
    Params.ClippingPlanes   = [-20, 20];
    Params.Perspective      = 1;
    Params.LightingOn       = 1;
    Params.Ambient          = [1 1 1 1];
    Params.Diffuse          = [1 1 1 1];
    Params.Specular         = [1 1 1 1];
    Params.Background       = [0.5 0.5 0.5 1];
end


%========================= OPEN GUI WINDOW ================================
Fig.Handle = figure;%(typecast(uint8('ENav'),'uint32'));               	% Assign GUI arbitrary integer        
if strcmp('ElectroNavInit', get(Fig.Handle, 'Tag')), return; end        % If figure already exists, return
Fig.FontSize = 12;
Fig.TitleFontSize = 14;
Fig.PanelYdim = 130;
Fig.Rect = [0 200 380 800];                                             % Specify figure window rectangle
set(Fig.Handle,     'Name','ElectroNav: 3D settings',...                % Open a figure window with specified title
                    'Tag','ElectroNav3D',...                            % Set figure tag
                    'Renderer','OpenGL',...                             % Use OpenGL renderer
                    'OuterPosition', Fig.Rect,...                       % position figure window
                    'NumberTitle','off',...                             % Remove figure number from title
                    'Resize', 'off',...                                 % Prevent resizing of GUI window
                    'Menu','none',...                                   % Turn off memu
                    'Toolbar','none');                                  % Turn off toolbars to save space
Fig.Background = get(Fig.Handle, 'Color');                              % Get default figure background color
Fig.FontSize = 10;                                                      % Set UI font size
Fig.Margin = 20;                                                        % Set margin between UI panels (pixels)                                 
Fig.Fields = fieldnames(Params);                                        % Get parameter field names


BoxPos{1} = [Fig.Margin,Fig.Rect(4)-220-Fig.Margin*2,Fig.Rect(3)-Fig.Margin*2, 220];         	% Set group controls positions
BoxPos{2} = [Fig.Margin,BoxPos{1}(2)-160-Fig.Margin/2,Fig.Rect(3)/2-Fig.Margin*1.5, 160];
BoxPos{3} = [Fig.Margin*2+BoxPos{2}(3),BoxPos{2}(2),BoxPos{2}(3),BoxPos{2}(4)];
BoxPos{4} = [Fig.Margin,BoxPos{2}(2)-160-Fig.Margin/2,Fig.Rect(3)-Fig.Margin*2, 160];
BoxPos{5} = [Fig.Margin,BoxPos{4}(2)-120-Fig.Margin/2,Fig.Rect(3)-Fig.Margin*2, 160];


% Logo= imread(fullfile('Documentation','ElectroNav_5.png'),'BackgroundColor',Fig.Background);
% LogoAx = axes('box','off','units','pixels','position', [120, 520, 260, 42],'color',Fig.Background);
% image(Logo);
% axis off

%=========================== SYSTEM PANEL =================================
Fig.SystemHandle = uipanel( 'Title','System Profile',...
                'FontSize',Fig.FontSize+2,...
                'BackgroundColor',Fig.Background,...
                'Units','pixels',...
                'Position',BoxPos{5},...
                'Parent',Fig.Handle); 
OS = computer;
MatlabVersion = version;
if exist('PsychtoolboxVersion', 'file')==2
    PTBversion = PsychtoolboxVersion;
else
    PTBversion = 'Not detected!';
end
OpenGL = opengl('data');

Ypos = BoxPos{5}(4)-Fig.Margin*2.5;
SystemValues = {OS, MatlabVersion,OpenGL.Renderer, OpenGL.Version, PTBversion};
SystemLabels = {'Operating system:','MATLAB version:','Graphics board:','OpenGL version:','PsychToolbox version:'};
for n = 1:numel(SystemLabels)
    h(n) = uicontrol('Style', 'text','String',SystemLabels{n},'Position', [Fig.Margin,Ypos,180,20],'Parent',Fig.SystemHandle,'HorizontalAlignment', 'left');
   	h(n+numel(SystemLabels)) = uicontrol('Style', 'text','String',SystemValues{n},'Position', [150,Ypos,180,20],'Parent',Fig.SystemHandle,'HorizontalAlignment', 'left');
 	Ypos = Ypos-25;
end
set(h, 'BackgroundColor',Fig.Background);



%% ======================== VIEWPORT GEOMETRY PANEL =======================
Fig.GeometryHandle = uipanel( 'Title','Scene Geometry',...
                'FontSize',Fig.FontSize+2,...
                'BackgroundColor',Fig.Background,...
                'Units','pixels',...
                'Position',BoxPos{1},...
                'Parent',Fig.Handle);  
            
Ypos = 170;           
%=================== STEREO MODE SELECTION
TipStr = 'Select a stereoscopic viewing method.';
RenderMethods = {'Monocular','Shutter glasses','Split-screen: top-bottom','Split-screen: bottom-top','Split-screen: free fusion','Split-screen: cross fusion',...
                 'Anaglyph: red-green','Anaglyph: green-red','Anaglyph: red-blue','Anaglyph: blue-red','Free fusion OSX',...
                 'PTB shutter glasses','Interleaved line','Interleaved column','Compressed HDMI'};
uicontrol(  'Style', 'text',...
            'String','Stereoscopic method:',...
            'Background',Fig.Background, ...
            'Position', [Fig.Margin,Ypos,150,20],...
            'parent',Fig.GeometryHandle,....
            'HorizontalAlignment', 'left');
uicontrol(  'Style', 'popup',...
            'Background', 'white',...
            'Tag', 'Method', ...
            'String', RenderMethods,...
            'Position', [178,Ypos,140,20],...
            'parent',Fig.GeometryHandle,...
            'TooltipString', TipStr,...
            'Callback', {@SetGeometry, 1});
Ypos = Ypos-25;
if strcmpi(PTBversion,'Not detected!')
    set(findobj('Tag','Method'), 'Value',1,'Enable','off');
else
    set(findobj('Tag','Method'), 'Value',Params.Stereomode+1);
end

ViewPortStrings = {'Viewing distance (cm):','Interpupillary distance (cm)','Screen dimensions (cm)','Frustum clipping planes (cm)'};
TipStr = {  'Set the viewing distance in centimetres (distance from observer to screen).',...
            'Set the observer''s interpupillary distance (distance between the eyes) in centimetres.',...
            'Set the physical dimensions of the display screen in centimetres (width x height)'...
            'Set the depth limits (near and far) beyond which objects will not be rendered'};
Tags = {'VD','IPD','ScreenDim','Zlims'};
DefaultAns = {Params.ViewingDist,Params.IPD,Params.ScreenDims,Params.ClippingPlanes};

%================== 
Indx = 2;
for n = 1:numel(ViewPortStrings)
    LabelH(n) = uicontrol(  'Style', 'text','String',ViewPortStrings{n},'Position', [Fig.Margin,Ypos,160,20],...
                'TooltipString', TipStr{n},'Parent',Fig.GeometryHandle,'HorizontalAlignment', 'left');
    if ismember(n,[3,4])
        Xpos = [180, 250];
        for i = 1:2
            uicontrol(  'Style', 'edit','Tag', Tags{n},'String', num2str(DefaultAns{n}(i)),'Position', [Xpos(i),Ypos,50,22],...
                        'TooltipString', TipStr{n},'Parent',Fig.GeometryHandle,'Callback', {@SetGeometry, Indx});
                    Indx = Indx+1;
        end
    else
        uicontrol(  'Style', 'edit','Tag', Tags{n},'String', num2str(DefaultAns{n}),'Position', [180,Ypos,120,22],...
                    'TooltipString', TipStr{n},'Parent',Fig.GeometryHandle,'Callback', {@SetGeometry, Indx});
                    Indx = Indx+1;
    end
    Ypos = Ypos-25;
end
set(LabelH, 'BackgroundColor',Fig.Background);


%================== PROJECTION
TipStr = 'Select the camera projection method.';
uicontrol(  'Style', 'text',...
            'String','Projection:',...
            'Background',Fig.Background, ...
            'Position', [Fig.Margin,Ypos,160,20],...
            'TooltipString', TipStr,...
            'Parent',Fig.GeometryHandle,...
            'HorizontalAlignment', 'left');
ProjHandle = uibuttongroup('visible','on','Background',Fig.Background,'units','pixels','Position',[180,Ypos-30,120,50],'parent', Fig.GeometryHandle);%, 'SelectionChangeFcn', {@SetGeometry, Indx});
ProjLabels = {'Perspective','Orthographic'};
for n = 1:numel(ProjLabels)
    uicontrol(  'Style', 'radio','Tag', 'Proj','String', ProjLabels{n},'Position', [10,5+(n-1)*20,120,20], 'TooltipString', TipStr,'Parent',ProjHandle,'Callback', {@SetGeometry, Indx});
    Indx = Indx+1;
    Ypos = Ypos-30;
end
% set(ProjHandle,'SelectedObject',Params.Perspective+1);

%================= LIGHTING PANEL
Fig.LightingHandle = uipanel( 'Title','Lighting',...
                'FontSize',Fig.FontSize+2,...
                'BackgroundColor',Fig.Background,...
                'Units','pixels',...
                'Position',BoxPos{2},...
                'Parent',Fig.Handle);  
Ypos = BoxPos{2}(4)-Fig.Margin*2-10;           
TipStr = 'Select light reflectance component colors.';
uicontrol(  'Style', 'togglebutton','Value',1,'String','Lighting on/off','Position', [Fig.Margin,Ypos,120,20],'Background',Fig.Background,...
            'TooltipString', TipStr,'Parent',Fig.LightingHandle,'HorizontalAlignment', 'left','Callback', {@SetLighting, 1});
LightingLabels = {'Ambient','Diffuse','Specular','Background'};
LightingDefaults = [Params.Ambient; Params.Diffuse; Params.Specular; Params.Background];     


% ShadeHandle = uibuttongroup('visible','off','Position',[180,Ypos,120,25],'parent', Fig.GeometryHandle, 'SelectionChangeFcn', {@FileCheck, 1});
Indx = 2;
for n = 1:numel(LightingLabels)
  	Ypos = Ypos-25;
    uicontrol(  'Style', 'radio',...
                'Background', Fig.Background,...
                'String', LightingLabels{n},...
                'Position', [15,Ypos,80,25],...
                'TooltipString', TipStr,...
                'Parent',Fig.LightingHandle,...
                'Callback', {@SetLighting, Indx});
    uicontrol(  'Style', 'pushbutton',...
                'Background', LightingDefaults(n,1:3),...
                'Tag', LightingLabels{n}, ...
                'String', '',...
                'Position',[110,Ypos,20,20],...
                'Parent',Fig.LightingHandle,...
                'Callback', {@SetLighting, Indx+1});
	Indx = Indx+2;
end



%================= MATERIALS PANEL
Fig.MaterialsHandle = uipanel( 'Title','Materials',...
                'FontSize',Fig.FontSize+2,...
                'BackgroundColor',Fig.Background,...
                'Units','pixels',...
                'Position',BoxPos{3},...
                'Parent',Fig.Handle);  
Ypos = BoxPos{3}(4)-Fig.Margin*2-10;           
TipStr = 'Select light reflectance component colors.';
uicontrol(  'Style', 'togglebutton','Value',1,'String','Materials on/off','Position', [Fig.Margin,Ypos,120,20],'Background',Fig.Background,...
            'TooltipString', TipStr,'Parent',Fig.MaterialsHandle,'HorizontalAlignment', 'left');
        
        
MaterialsLabels = {'Ambient','Diffuse','Specular','Background'};
        
% ShadeHandle = uibuttongroup('visible','off','Position',[180,Ypos,120,25],'parent', Fig.GeometryHandle, 'SelectionChangeFcn', {@FileCheck, 1});
for n = 1:numel(MaterialsLabels)
  	Ypos = Ypos-25;
    uicontrol(  'Style', 'radio',...
                'Background', Fig.Background,...
                'String', MaterialsLabels{n},...
                'Position', [15,Ypos,80,25],...
                'TooltipString', TipStr,...
                'Parent',Fig.MaterialsHandle);
    uicontrol(  'Style', 'pushbutton',...
                'Background', [0.3 0.3 0.3],...
                'Tag', MaterialsLabels{n}, ...
                'String', '',...
                'Position',[110,Ypos,20,20],...
                'Parent',Fig.MaterialsHandle,...
                'Callback', {@SetColor, n});
end
    

%================= OPTIONS PANEL
uicontrol(  'Style', 'pushbutton',...
            'String','Load',...
            'parent', Fig.Handle,...
            'tag','Load',...
            'units','pixels',...
            'Position', [Fig.Margin,20,100,30],...
            'TooltipString', 'Use current inputs',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 1});   
uicontrol(  'Style', 'pushbutton',...
            'String','Save',...
            'parent', Fig.Handle,...
            'tag','Save',...
            'units','pixels',...
            'Position', [140,20,100,30],...
            'TooltipString', 'Save current inputs to file',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 2});    
uicontrol(  'Style', 'pushbutton',...
            'String','Continue',...
            'parent', Fig.Handle,...
            'tag','Continue',...
            'units','pixels',...
            'Position', [260,20,100,30],...
            'TooltipString', 'Exit',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 3});         

hs = guihandles(Fig.Handle);                                % get UI handles
guidata(Fig.Handle, hs);                                    % store handles
set(Fig.Handle, 'HandleVisibility', 'callback');            % protect from command line
drawnow;
% uiwait(Fig.Handle);
ParamsOut = Params;




%% ========================= UICALLBACK FUNCTIONS =========================
    function SetGeometry(hObj, Evnt, Indx)
        switch Indx
            case 1          %============== Stereomode
                Params.Stereomode = get(hObj,'Value')-1;
            case 2
                Params.ViewingDist = str2num(get(hObj,'String'));
            case 3
                Params.IPD = str2num(get(hObj,'String'));
            case 4
                Params.ScreenDim(1) = str2num(get(hObj,'String'));
            case 5
                Params.ScreenDim(2) = str2num(get(hObj,'String'));
            case 6
                Params.ClippingPlanes(1) = str2num(get(hObj,'String'));
                if Params.ClippingPlanes(1)>0
                    Params.ClippingPlanes(1) = -Params.ClippingPlanes(1);
                end
            case 7
                Params.ClippingPlanes(2) = str2num(get(hObj,'String'));
                if Params.ClippingPlanes(2)<0
                    Params.ClippingPlanes(2) = -Params.ClippingPlanes(2);
                end
            case 8
                Params.Perspective = 0;
            case 9
                Params.Perspective = 1;
        end
        Params

    end


    function SetLighting(hObj, Evnt, Indx)
        switch Indx
            
            case 1
                Params.LightingOn = get(hObj,'Value');
            case {2,4,6,8}
                
            case {3,5,7,9}
                Color = uisetcolor;
                set(hObj, 'Background', Color);
                Params.Ambient = [Color, 1];
        end
    end


    %==================== OPTIONS
    function OptionSelect(Obj, Event, Indx)
        switch Indx
            case 1      %================ Continue
                RequiredFieldsIndx = [3,4,6];
                RequiredFieldsAbsent = zeros(1,numel(RequiredFieldsIndx));
                for r = 1:numel(RequiredFieldsIndx)
                    RequiredFields{r} = eval(sprintf('Params.%s;',Fig.Fields{RequiredFieldsIndx(r)}));
                    if isempty(RequiredFields{r})
                        RequiredFieldsAbsent(r) = 1;
                    end
                end
                if any(RequiredFieldsAbsent)
                    h = errordlg('Required inputs are missing!', 'Insufficient inputs!','modal');
                    uiwait(h);
                    return;
                else
                    close(Fig.Handle);      % Close GUI figure
                    return;
                end

            case 2      %================ Save
                if strcmp(Params.SubjectID,'')
                    return;
                end
                Indx = structfind(Params.Defaults,'SubjectID',Params.SubjectID);
                for f = 1:numel(Fig.Fields)
                    eval(sprintf('Params.Defaults(Indx).%s = Params.%s;',Fig.Fields{f},Fig.Fields{f}));
                end
                Defaults = Params.Defaults;
                [Filename, Pathname, Indx] = uiputfile('*.mat','Save default parameters file','DefaultParamsFile.mat');
                if Filename == 0
                    return;
                end
                save(fullfile(Pathname, Filename),'Defaults');
                msgbox('Default parameters saved!','Saved');

            case 3      %================ Cancel
                ParamsOut = [];         % Clear params
                close(Fig.Handle);      % Close GUI figure
                return;
        end
    end

end