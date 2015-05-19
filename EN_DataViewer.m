function [Contact] = EN_DataViewer(DataPath)

%=========================== EN_DataViewer.m ==============================
% This function loads analysed physiology data from multiple prior sessions 
% and plots the results in 3D space, overlaid with segmented neural atlas
% stuctures for visualization of spatial patterns.
%
% INPUTS:
%   The input 'DataPath' should be the directory of .mat files containing 
%   the following structure:
%
%       Contact.Dates:      A 1 x d cell array of date strings of recording
%                           session(s), in any standard Matlab date format.
% 
%       Contact.DateFormat: Format off date strings provided e.g.'yyyymmdd';
% 
%       Contact.XYZ:       	a [Session# x Channel# x XYZ] matrix
% 
%       Contact.CellIndxData: an N x 4 matrix with columns containing: 
%                           [CellID#, Session#, Channel#, Cell#]
% 
%       Contact.ColorVals: 	A d x n matrix containing a scalar value to specify the 
%                           color for each data point via a colormap.
%       Contact.rad:        A d x n matrix containing a scalar value to specify the
%                           radius of each data point in millimetres.
%       Contact.Alpha:      A d x n matrix containing a scalar value to specify the
%                           opacity (1-transparency) of each data point.
%   
% REVISIONS:
%   27/04/2014 - Written by APM
%   19/03/2015 - GUI developed
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, � Copyleft 2015, GNU General Public License
%==========================================================================

clear all;
[root, temp] = fileparts(mfilename('fullpath'));
addpath(genpath(root));
global Fig Contact Structures Data

%==================== LOAD PHYSIOLOGY RESULTS DATA
if nargin == 0
    Data.Dir = uigetdir(root, 'Select data directory');
else
    Data.Dir = DataPath;
end                         
Data.Files = wildcardsearch(Data.Dir,'*.mat');  
for d = 1:numel(Data.Files)
    [a,Data.Filenames{d},c] = fileparts(Data.Files{d});
end
[Data.Selected,ok] = listdlg('PromptString','Select data', 'SelectionMode','single','ListString',Data.Filenames,'ListSize',[300, 250]);
load(Data.Files{Data.Selected});                                                          	% Load selected data file
if ~exist('Contact','var')                                                                  % Check that file contained Contact structure
	error(sprintf('No ''contact'' structure was found in %s!', Data.Files{Data.Selected})); 
end
if ~any(isfield(Contact,{'Dates','ColorVals','Alpha','rad','XYZ','CellIndxData'}))          % Check that structure fontains required fields
    error('Insufficient data fields supplied!');
end

%========================= OPEN FIGURE WINDOW =============================
LoadingFig = EN_About(1);                                                       % Open loading screen
Fig.scnsize = get(0,'ScreenSize');                                              % Get screen resolution
Fig.Rect = [0 0 Fig.scnsize(3), Fig.scnsize(4)];                                % Set figure winow to fullscreen
Fig.FontSize = 14;                                                              % Set defualt font size
Fig.Background = repmat(0.75,[1,3]);                                            % Set figure background color  
Fig.InputBackground = repmat(0.85,[1,3]);                                       % Set input box background color
Fig.AxesBkgColor = repmat(0.75,[1,3]);                                          % Set axes background color
Fig.Handles.Figure = figure('Name',sprintf('ElectroNav%c - Data Viewer',char(169)),... 	% Open a figure window with specified title
                    'Color',Fig.Background,...                                  % Set the figure window background color
                    'Renderer','OpenGL',...                                     % Use OpenGL renderer
                    'Position', Fig.Rect,...                                    % position figure window to fit fullscreen
                    'visible','off',...                                         % Figure remains invisible until complete
                    'menu','none','toolbar','none',...                          % Remove toolbar etc
                    'NumberTitle','off',...                                     % Remove figure number from title
                    'IntegerHandle','off');                                     % Don't use integer handles
Fig.Colormap = 'parula';                                                        % Set default colormap
Fig.HighlightColor = [1 0 0];                                                  	% Set color for highlighted contacts
Fig.HighlightAlpha = 0.3;                                                       % Set alpha transparency for outer highlight volume
Fig.HighlightRad = 0.5;                                                         % Set radius (mm) of outer highlight volume
Fig.View = [220, 30];                                                           % Set default camera angle [azimuth, elevation] (degrees)
Fig.SphereRadii = [0.05, 0.15];                                               	% Set the range of radii (mm) of contacts
Fig.SphereN = 20;                                                             	% Set number of segments per sphere (redcude for faster rendering, increase for higher quality)
Fig.Xlim = [-16 0];                                                             % Set default axis limits for 3D plot
Fig.Ylim = [-20,-8];
Fig.Zlim = [-8,8];
Fig.MarkerSize = 10;                                                            % Default marker size for 2D plots
Fig.MarkerColor = [0 0 1];                                                      % Default marker color for 2D plots
Fig.Data.InvertAlpha = 0;                                                       % Default alpha transparency is NOT inverted (higher alpha = more opaque)


%================ PLOT 3D DATA
RadRange = [min(Contact.rad), max(Contact.rad)];                                    % Get range of radius values
if numel(unique(Contact.rad))>1
    Contact.NormRad = (Contact.rad-RadRange(1))/ diff(RadRange);                  	% Normalize radius values (0-1)
    Contact.NormRad = Fig.SphereRadii(1)+(Contact.NormRad*diff(Fig.SphereRadii));  	% Scale radius values to specified range
else
    Contact.NormRad = ones(size(Contact.rad))*Fig.SphereRadii(2);                   % If all radii supplied are same, use max
end
AlphaRange = [min(Contact.Alpha), max(Contact.Alpha)];                              % Get range of alpha values
if diff(AlphaRange)==0
    Contact.NormAlpha = ones(size(Contact.Alpha));                                  
else
    Contact.NormAlpha = (Contact.Alpha-AlphaRange(1))/ diff(AlphaRange);         	% Normalize radius values (0-1)
end
Fig.Handles.MainAx = axes('units','normalized','position',[0.05,0.1,0.5,0.85]);
for d = [unique(Contact.CellIndxData(:,2))]'                                      	% For each session date...
    SessionIndx = find(Contact.CellIndxData(:,2)==d);                               % Get indices of current session
    for ch = [unique(Contact.CellIndxData(SessionIndx,3))]'                     	% For each contact location...
        ChannelIndx = SessionIndx(find(Contact.CellIndxData(SessionIndx,3)==ch));   % Get indices of current channel
        NoCells = numel(unique(Contact.CellIndxData(ChannelIndx,4)));              	% Find the number of cells
        for c = unique(Contact.CellIndxData(ChannelIndx,4))'                     	% For each cell...
            CellIndx = ChannelIndx(find(Contact.CellIndxData(ChannelIndx,4)==c));   % Get the cell index
            if NoCells==1                                                        	% If only 1 cell was recorded on this channel....
                [X,Y,Z] = ellipsoid(Contact.XYZ(d,ch,1),Contact.XYZ(d,ch,2),Contact.XYZ(d,ch,3),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Fig.SphereN);
            elseif NoCells>1
                theta = (360/NoCells)*c;
                Offset = Fig.SphereRadii(2)*[cosd(theta), sind(theta)];
                [X,Y,Z] = ellipsoid(Contact.XYZ(d,ch,1)+Offset(1),Contact.XYZ(d,ch,2)+Offset(2),Contact.XYZ(d,ch,3),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Fig.SphereN);
            end
            Contact.Handles(CellIndx) = surface(X,Y,Z,repmat(Contact.ColorVals(CellIndx),size(Z)),'EdgeColor','none','CDataMapping','scaled','FaceAlpha',Contact.NormAlpha(CellIndx));
            hold on;
        end
    end
end
set(Contact.Handles,'ButtonDownFcn',@ContactSelect);                                % Set callback function for mouse selection of 3D data points    

%================ Set plot appearance
shading interp; 
lighting phong; 
view(Fig.View(1), Fig.View(2));
Fig.Handles.Light(1) = camlight('left'); 
Fig.Handles.Light(2) = camlight('headlight'); 
grid on;
axis equal tight vis3d;
xlabel('Medial-Lateral (mm)','fontsize',18);
ylabel('Anterior-Posterior (mm)','fontsize',18);
zlabel('Inferior-Superior (mm)','fontsize',18);
set(gca,'fontsize',Fig.FontSize,'color',Fig.AxesBkgColor);                      
set(gca,'xlim',Fig.Xlim','ylim',Fig.Ylim,'zlim',Fig.Zlim);                      % Set axes limits to defaults
colormap(Fig.Colormap);                                                         % Set colormap to requested default
Fig.cbh = colorbar('EastOutside','units','pixels','position',[80 200 25 500]);  % Plot colorbar scale
Fig.cbh.Label.String = Contact.DataVariable;                                    % Set axis label on colorbar
Fig.cbh.Label.FontSize = 18;                                                    % Set fontsize for axis label
Fig.cbh.FontSize = 12;                                                          
Fig.cbh.ButtonDownFcn = @SelectColormap;                                        % Set callback function for colorbar (select new colormap)
Contact.Indx = 1;                                                               % Currently selected contact/cell defaults to 1
Fig.ZoomH = zoom;                                                               % Create zoom object

%================ Plot structures
Structures = Plot3DMeshes;
Structures.Materials.Ambient = 0.3;                              
Structures.Materials.Diffuse = 0.5;          
Structures.Materials.Specular = 0.4;         
Structures.Materials.SpecExp = 6;            
Structures.Materials.SpecCol = 1; 
material(cell2mat(struct2cell(Structures.Materials))');

Structures.CurrentStructure = 1;
Structures.Smoothing = zeros(1, numel(Structures.Names));
Structures.On = zeros(1, numel(Structures.Names));
Structures.On([2,5,6,8,9]) = 1;                             
Structures.Wire = ones(1, numel(Structures.Names));
for s = 1:numel(Structures.Handles)
    if Structures.On(s)==1
        set(Structures.Handles{s},'visible','on');
    else
        set(Structures.Handles{s},'visible','off');
    end
end


%% =========================== ADD GUI PANELS ===============================
Fig.Handles.OuterPannel = uipanel('BackgroundColor',Fig.Background,'Units','normalized','Position',[0.58,0.05,0.4,0.9]);

PanelNames = {'Selected cell', 'Atlas structures', 'Options', 'Data'};
BoxPos(1,:) = [400, 600, 260, 240];
BoxPos(2,:) = [BoxPos(1,1), BoxPos(1,2)-BoxPos(1,4)-20, BoxPos(1,3), 240];
BoxPos(3,:) = [BoxPos(1,1), BoxPos(2,2)-BoxPos(2,4)-20, BoxPos(1,3), 240];
BoxPos(4,:) = [20, BoxPos(3,2), 360, sum(BoxPos(1:3,4))+40];
Fig.Handles.UIpannel(1) = uipanel('Title',PanelNames{1},'FontSize',14,'BackgroundColor',Fig.Background,'Units','pixels','Position',BoxPos(1,:),'Parent',Fig.Handles.OuterPannel);
Fig.Handles.UIpannel(2) = uipanel('Title',PanelNames{2},'FontSize',14,'BackgroundColor',Fig.Background,'Units','pixels','Position',BoxPos(2,:),'Parent',Fig.Handles.OuterPannel);
Fig.Handles.UIpannel(3) = uipanel('Title',PanelNames{3},'FontSize',14,'BackgroundColor',Fig.Background,'Units','pixels','Position',BoxPos(3,:),'Parent',Fig.Handles.OuterPannel);
Fig.Handles.UIpannel(4) = uipanel('Title',PanelNames{4},'FontSize',14,'BackgroundColor',Fig.Background,'Units','pixels','Position',BoxPos(4,:),'Parent',Fig.Handles.OuterPannel);


%========================= CELL INFO PANEL
Fig.Labels = {'Date','[x,y,z] (mm)','Contact #','Cell #','Color value','Alpha value','Radius value'};
Fig.Tags = {'Date','x','y','z','ContactNo','CellNo','Color', 'Alpha','Radius'};
Fig.Formats = {'%s','%.1f','%.1f','%.1f','%d','%d','%.2f','%.2f','%.2f'};
Ypos = BoxPos(1,4)-20;
f = 1;
for n = 1:numel(Fig.Labels)
    Ypos = Ypos-30;
    Fig.Handles.CellLabel(n) = uicontrol('Style','Text','String',Fig.Labels{n},'HorizontalAlignment','Left','pos',[10, Ypos, 100, 25],'parent',Fig.Handles.UIpannel(1));
 	if n == 1
        Fig.Handles.CellInput(f) = uicontrol('Style','popup','String',Contact.Dates,'tag',Fig.Tags{f},'pos',[90, Ypos, 150, 25],'Callback',{@ContactInput,f},'parent',Fig.Handles.UIpannel(1));      
        f = f+1;
 	elseif n == 2
        Xpos = 90;
        for i = 1:3
           Fig.Handles.CellInput(f) = uicontrol('Style','Edit','String','','tag',Fig.Tags{f},'pos',[Xpos, Ypos, 48, 25],'Callback',{@ContactInput,f},'parent',Fig.Handles.UIpannel(1)); 
            Xpos = Xpos+50;
            f = f+1;
        end
    else
        Fig.Handles.CellInput(f) = uicontrol('Style','Edit','String','','tag',Fig.Tags{f},'pos',[90, Ypos, 150, 25],'Callback',{@ContactInput,f},'parent',Fig.Handles.UIpannel(1));
    	f = f+1;
    end
end
UpdateUIFields(1);      % Set default information to first cell


%================ STRUCTURE OVERLAY PANEL
Fig.Struct.ButtonDim = [100 25];
Fig.Struct.LabelStrings = {'Current structure','Color','Opacity','Smoothing','',''};
Fig.Struct.InputType = {'popupmenu','PushButton','slider','slider', 'checkbox','checkbox'};
Fig.Struct.InputStrings = {Structures.Names, [], [], [], 'Visible','Wireframe'};
Fig.Struct.InputValue = {Structures.CurrentStructure, [], Structures.Opacity(Structures.CurrentStructure), Structures.Smoothing(Structures.CurrentStructure), Structures.On(Structures.CurrentStructure), Structures.Wire(Structures.CurrentStructure)};
Fig.Struct.ButtonPos = [repmat(10,[numel(Fig.Struct.LabelStrings),1]), [0:30:((numel(Fig.Struct.LabelStrings)-1)*30)]'+30];
Fig.Struct.ButtonPos = Fig.Struct.ButtonPos(end:-1:1,:);
for i = 1:numel(Fig.Struct.LabelStrings)
    Fig.Handles.StructLabel(i) = uicontrol('Style','text', 'string', Fig.Struct.LabelStrings{i},'HorizontalAlignment','Left', 'pos', [Fig.Struct.ButtonPos(i,:), Fig.Struct.ButtonDim],'parent',Fig.Handles.UIpannel(2));
    Fig.Handles.StructInput(i) = uicontrol('Style',Fig.Struct.InputType{i},'String',Fig.Struct.InputStrings{i},'value',Fig.Struct.InputValue{i}, 'pos',[Fig.Struct.ButtonPos(i,:)+[Fig.Struct.ButtonDim(1), 0], Fig.Struct.ButtonDim+[20 0]],'parent',Fig.Handles.UIpannel(2),'Callback',{@StructView,i}); 
end
set(Fig.Handles.StructLabel(3),'String',sprintf('Opacity = %d %%', Structures.Opacity(Structures.CurrentStructure)*100));
% set(Fig.Handles.StructLabel(2),'String',sprintf('Smoothing = %.0f mm', Fig.Struct.sigma));
set(Fig.Handles.StructInput(2), 'BackgroundColor', Structures.Colors(Structures.CurrentStructure, :));


%======================== USER ACTIONS PANEL
Fig.Labels = {'Manually rotate','Camera position','Light position','Set axes limits','Save settings','Load settings','Zoom','Export figure','Stereo render'};
Fig.ButtonSyles = {'Togglebutton','Pushbutton','Pushbutton','Pushbutton','Pushbutton','Pushbutton','Togglebutton','Pushbutton','Pushbutton'};
Ypos = (BoxPos(3,4)-50)+ (0:-30:(-30*5));
for n = 1:numel(Fig.Labels)
    if n <=6 
        ButtonPos = [20, Ypos(n), 100, 25];
    else
        ButtonPos = [140, Ypos(n-6), 100, 25];
    end
    Fig.Handles.OptionsInput(n) = uicontrol('Style',Fig.ButtonSyles{n},...
                                'String',Fig.Labels{n},...
                             	'HorizontalAlignment','Left',...
                              	'pos',ButtonPos,...
                                'Callback',{@OptionsInput,n},...
                              	'parent',Fig.Handles.UIpannel(3));
end
set(Fig.Handles.OptionsInput([1,7]),'value',0);                                 % Zoom and rotate default to 'off'

%======================== DATA PANEL
Fig.Data.AxisLabels = {'Medial-lateral','Anterior-posterior','Ventral-Dorsal'};
Fig.AxisSelected = 3;
Fig.Data.LabelStrings = {'Filename:','Number of sessions:','Number of neurons:','Selected axis:','Color thresholds:','Alpha thresholds:','Invert alpha:'};
Fig.Data.InputType = {'popupmenu','text','text','popupmenu','edit','edit','checkbox','edit','edit'};
ColorThresh = [min(Contact.ColorVals), max(Contact.ColorVals)];
AlphaThresh = [min(Contact.Alpha), max(Contact.Alpha)];
Fig.Data.InputStrings = {Data.Filenames, num2str(numel(Contact.Dates)),num2str(size(Contact.CellIndxData,1)),Fig.Data.AxisLabels, num2str(ColorThresh(1)), num2str(AlphaThresh(1)), 'Inverted',num2str(ColorThresh(2)), num2str(AlphaThresh(2))};
Fig.Data.InputValue = {Data.Selected, [], [], Fig.AxisSelected, [], [], Fig.Data.InvertAlpha, [], []};
Ypos = (0:-30:(-30*(numel(Fig.Data.LabelStrings)+1))) + BoxPos(4,4)-50;
Ypos([8,9]) = Ypos([5,6]);
InputXpos = [140 140 140 140 140 140 140 240 240];
InputWidth = [180 180 180 180 80 80 180 80 80];
for i = 1:numel(Fig.Data.LabelStrings)+2
    if i<=7
        Fig.Handles.DataLabel(i) = uicontrol('Style','text', 'string', Fig.Data.LabelStrings{i},'HorizontalAlignment','Left', 'pos', [10, Ypos(i), 120, 25],'parent',Fig.Handles.UIpannel(4));
    end
	Fig.Handles.DataInput(i) = uicontrol('Style',Fig.Data.InputType{i},'String',Fig.Data.InputStrings{i},'value',Fig.Data.InputValue{i}, 'pos',[InputXpos(i), Ypos(i), InputWidth(i), 25],'parent',Fig.Handles.UIpannel(4),'Callback',{@DataView,i});
end
set(Fig.Handles.DataInput(2:3), 'HorizontalAlignment','Left', 'Callback', []);
set(Fig.Handles.DataInput(7), 'value', Fig.Data.InvertAlpha);

%========= SET FIGURE/ PANNEL COLORS
Fig.Handles.AllPannels = [Fig.Handles.UIpannel, Fig.Handles.OuterPannel, Fig.Handles.CellLabel, Fig.Handles.OptionsInput, Fig.Handles.DataLabel, Fig.Handles.DataInput, Fig.Handles.StructLabel, Fig.Handles.StructInput([5,6])];
Fig.Handles.AllInputs = [Fig.Handles.CellInput, Fig.Handles.DataInput([5,6,8,9])];
set(Fig.Handles.AllPannels, 'BackgroundColor',Fig.Background);
set(Fig.Handles.AllInputs, 'BackgroundColor',Fig.InputBackground);
set(Fig.Handles.Figure, 'Color', Fig.Background);
set(Fig.Handles.MainAx, 'Color', Fig.InputBackground);


%======== PLOT DATA HISTOGRAM
Fig.Handles.DataAx(1) = axes('units','pixels','position',[60, 40, 270, 200],'Parent',Fig.Handles.UIpannel(4));
Fig.Data.Hist = histogram(Contact.ColorVals, round(numel(Contact.ColorVals)/2));
Fig.Data.Hist.EdgeColor = 'none';
hold on;
axis tight
grid on;
ylabel('Number of cells');
xlabel('Data value');

Fig.Handles.DataAx(2) = axes('units','pixels','position',[60, 300, 270, 200],'Parent',Fig.Handles.UIpannel(4));
Allpos = reshape(Contact.XYZ(:,:,Fig.AxisSelected)',[size(Contact.XYZ,1)*size(Contact.XYZ,2), 1]);
% t= linspace(0, 2*pi, 20);
for i = 1:size(Contact.CellIndxData,1)
    PosIndx(i) = ((Contact.CellIndxData(i,2)-1)*size(Contact.XYZ,2))+Contact.CellIndxData(i,3);
    Fig.Data.Spatial.ph(i) = plot(Allpos(PosIndx(i)), Contact.ColorVals(i),'.b','ButtonDownFcn',{@DataSelect, i});
%     Fig.Data.Spatial.ph(i) = patch(sin(t)+Allpos(PosIndx(i)), cos(t)+Contact.ColorVals(i), 'b','edgecolor','none','facecolor', Fig.MarkerColor,'FaceAlpha', Contact.Alpha(i));
    hold on;
end
set(Fig.Data.Spatial.ph,'MarkerSize',Fig.MarkerSize,'MarkerFaceColor',Fig.MarkerColor,'MarkerEdgeColor',Fig.MarkerColor);
grid on;
ylabel('Data value');
xlabel([Fig.Data.AxisLabels{Fig.AxisSelected},' (mm)']);
set(Fig.Handles.DataAx,'tickdir','out','box','off','fontsize',12);


close(LoadingFig);                              % Close the 'Loading...' message window
set(Fig.Handles.Figure,'visible','on');         % Make main figure window visible
end


%% ====================== DISPLAY SELECTED CELL DATA ======================
function ContactSelect(objectHandle, eventData)
global Contact Fig
    axesHandle  = get(objectHandle,'Parent');                               % Get axes handle
    Contact.Indx = find(Contact.Handles==objectHandle);                     % Get cell index number
    Contact.Current.DateIndx = Contact.CellIndxData(Contact.Indx,2);
    Contact.Current.ChannelIndx = Contact.CellIndxData(Contact.Indx,3);  
    Contact.Current.CellIndx = Contact.CellIndxData(Contact.Indx,4);
    UpdateHighlight(Contact);
    UpdateUIFields(Contact.Indx);
end

%==================== UPDATE VALUES IN UI FIELDS ==========================
function UpdateUIFields(ContactIndx)
global Contact Fig
 	Contact.Current.DateIndx = Contact.CellIndxData(Contact.Indx,2);
    Contact.Current.ChannelIndx = Contact.CellIndxData(Contact.Indx,3);  
    Contact.Current.CellIndx = Contact.CellIndxData(Contact.Indx,4);
    Coordinates = Contact.XYZ(Contact.Current.DateIndx,Contact.Current.ChannelIndx,:);  	% Get coordinates of selected contact
    Inputs = {Contact.Current.DateIndx, Coordinates(1),Coordinates(2),Coordinates(3),...
              Contact.Current.ChannelIndx,Contact.Current.CellIndx, Contact.ColorVals(ContactIndx),Contact.Alpha(ContactIndx),Contact.rad(ContactIndx)};
    for n = 1:numel(Fig.Tags)                                       % For each GUI field tag...
        if n==1                                                     % For date popup menu...
            set(findobj('Tag',Fig.Tags{n}), 'Value', Inputs{n});
        else
            Inputs{n} = sprintf(Fig.Formats{n},Inputs{n});        	% Convert numbers to correct format string
            set(findobj('Tag',Fig.Tags{n}), 'String', Inputs{n});   % Set field to new value
        end
    end
end

%==================== UPDATE POSITION OF CONTACT HIGHLIGHT
function UpdateHighlight(Contact)
global Fig
    Coordinates = squeeze(Contact.XYZ(Contact.Current.DateIndx,Contact.Current.ChannelIndx,:)); 	% Get coordinates of selected contact
    set(Contact.Handles,'EdgeColor','none');                                            % Turn edges off for all contacts
    set(Contact.Handles(Contact.Indx),'EdgeColor',Fig.HighlightColor);                	% Turn edges of selected contact red
    if isfield(Fig,'Highlight') & ishandle(Fig.Highlight)
        delete(Fig.Highlight);                                                          % Delete 3D surface marker
    end
    set(Fig.Data.Spatial.ph, 'markerfacecolor', Fig.MarkerColor, 'markeredgecolor', Fig.MarkerColor, 'MarkerSize', Fig.MarkerSize);                             % Restore defaults for all markers in 2D plot
    set(Fig.Data.Spatial.ph(Contact.Indx),'markerfacecolor', Fig.HighlightColor, 'markeredgecolor', Fig.HighlightColor, 'MarkerSize', Fig.MarkerSize+20);       % Apply highlight to selected marker
    uistack(Fig.Data.Spatial.ph(Contact.Indx), 'top');                                                                                                          % Move selected datum to top layer
	axes(Fig.Handles.MainAx);                                                                                                                                   % Select the main 3D axes
    [X,Y,Z] = ellipsoid(Coordinates(1),Coordinates(2),Coordinates(3),Fig.HighlightRad,Fig.HighlightRad,Fig.HighlightRad,50);                                    % Get coordinates for a larger highlight sphere
    Fig.Highlight = surface(X,Y,Z,repmat(mean(Contact.ColorVals), size(Z)),'EdgeColor','none','FaceColor',Fig.HighlightColor,'FaceAlpha',Fig.HighlightAlpha);   % Plot the highlight sphere
    if isfield(Fig.Data,'DistHandle') & ishandle(Fig.Data.DistHandle)
        set(Fig.Data.DistHandle, 'xdata', repmat(Contact.ColorVals(Contact.Indx),[1,2]));
    else
        axes(Fig.Handles.DataAx(1));
        Fig.Data.DistHandle = plot(repmat(Contact.ColorVals(Contact.Indx),[1,2]),ylim, '-r', 'linewidth',2);
    end
end

%====================== DISPLAY SELECTED CELL DATA
function DataSelect(objectHandle, eventData, Indx)
global Contact Fig
    Contact.Indx = Indx;                                                        % Get cell index number
    axesHandle  = get(objectHandle,'Parent');                                   % Get axes handle
    Contact.Current.DateIndx = Contact.CellIndxData(Contact.Indx,2);
    Contact.Current.ChannelIndx = Contact.CellIndxData(Contact.Indx,3);  
    Contact.Current.CellIndx = Contact.CellIndxData(Contact.Indx,4);
    UpdateHighlight(Contact);
    UpdateUIFields(Contact.Indx);
end

%====================== DISPLAY SELECTED CELL DATA
function SelectColormap(objectHandle, eventData)
global Fig
    AllCMaps = {'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'};
    [SelectedCMap,ok] = listdlg('PromptString','Select a colormap', 'SelectionMode','single','ListString',AllCMaps,'ListSize',[250, 120]);
    if ok == 1
        axes(Fig.Handles.MainAx);
        Fig.Colormap = AllCMaps{SelectedCMap};  
        colormap(Fig.Colormap);
    end
end

%========================== LOAD NEW DATA =================================
function LoadNewData(Contact)
global Fig Contact
    %=========== Update GUI information
    set(Fig.Handles.CellInput(1),'string',Contact.Dates);                                 % List new session dates
    Contact.Indx = 1;                                                               % Reset contact index
    UpdateUIFields(Contact.Indx);                                                   % Update contact pannel info
  	set(Fig.Handles.DataInput(2),'string',num2str(numel(Contact.Dates)));          	% Update data pannel info
    set(Fig.Handles.DataInput(3),'string',num2str(size(Contact.CellIndxData,1)));    
    Thresh = [min(Contact.ColorVals), max(Contact.ColorVals)];                      % Get minimum and maximum data values
    set(Fig.Handles.DataInput(5),'string',num2str(Thresh(1)));                    	% Set lower threshold value
    set(Fig.Handles.DataInput(8),'string',num2str(Thresh(2)));                     	% Set upper threshold value
 	AlphaThresh = [min(Contact.Alpha), max(Contact.Alpha)];                         % Get minimum and maximum data values
    set(Fig.Handles.DataInput(6),'string',num2str(AlphaThresh(1)));              	% Set lower threshold value
    set(Fig.Handles.DataInput(9),'string',num2str(AlphaThresh(2)));               	% Set upper threshold value

    %=========== Plot new data to graphs
    axes(Fig.Handles.DataAx(1));                                                           % Select histogram axes
    Fig.Data.Hist = histogram(Contact.ColorVals, round(numel(Contact.ColorVals)/2));% Plot histogram
    Fig.Data.Hist.EdgeColor = 'none';
    axes(Fig.Handles.DataAx(2));                                                         	% Select data axes
    Allpos = reshape(Contact.XYZ(:,:,Fig.AxisSelected)',[size(Contact.XYZ,1)*size(Contact.XYZ,2), 1]);
    for i = 1:size(Contact.CellIndxData,1)
        PosIndx(i) = ((Contact.CellIndxData(i,2)-1)*size(Contact.XYZ,2))+Contact.CellIndxData(i,3);
    end
    Fig.Data.Spatial.ph = plot(Allpos(PosIndx), Contact.ColorVals,'.b');            % Plot position versus data
 	AllLims = [Fig.Xlim; Fig.Ylim; Fig.Zlim];
    set(Fig.Handles.DataAx(1), 'xlim', Thresh);
    set(Fig.Handles.DataAx(2), 'ylim', Thresh);
 	set(Fig.Handles.DataAx(2),'xlim', AllLims(Fig.AxisSelected,:));
    
    %=========== Normalize data
    RadRange = [min(Contact.rad), max(Contact.rad)];                                    % Get range of radius values
    if numel(unique(Contact.rad))>1
        Contact.NormRad = (Contact.rad-RadRange(1))/ diff(RadRange);                  	% Normalize radius values (0-1)
        Contact.NormRad = Fig.SphereRadii(1)+(Contact.NormRad*diff(Fig.SphereRadii));  	% Scale radius values to specified range
    else
        Contact.NormRad = ones(size(Contact.rad))*Fig.SphereRadii(2);                   % If all radii supplied are same, use max
    end
    AlphaRange = [min(Contact.Alpha), max(Contact.Alpha)];                              % Get range of alpha values
    if diff(AlphaRange)==0
        Contact.NormAlpha = ones(size(Contact.Alpha));                                  
    else
        Contact.NormAlpha = (Contact.Alpha-AlphaRange(1))/ diff(AlphaRange);         	% Normalize radius values (0-1)
    end
    
    %=========== Plot new 3D data points
    axes(Fig.Handles.MainAx);
    for d = [unique(Contact.CellIndxData(:,2))]'                                      	% For each session date...
        SessionIndx = find(Contact.CellIndxData(:,2)==d);                               % Get indices of current session
        for ch = [unique(Contact.CellIndxData(SessionIndx,3))]'                     	% For each contact location...
            ChannelIndx = SessionIndx(find(Contact.CellIndxData(SessionIndx,3)==ch));   % Get indices of current channel
            NoCells = numel(unique(Contact.CellIndxData(ChannelIndx,4)));              	% Find the number of cells
            for c = unique(Contact.CellIndxData(ChannelIndx,4))'                     	% For each cell...
                CellIndx = ChannelIndx(find(Contact.CellIndxData(ChannelIndx,4)==c));   % Get the cell index
                if NoCells==1                                                        	% If only 1 cell was recorded on this channel....
                    [X,Y,Z] = ellipsoid(Contact.XYZ(d,ch,1),Contact.XYZ(d,ch,2),Contact.XYZ(d,ch,3),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Fig.SphereN);
                elseif NoCells>1
                    theta = (360/NoCells)*c;
                    Offset = Fig.SphereRadii(2)*[cosd(theta), sind(theta)];
                    [X,Y,Z] = ellipsoid(Contact.XYZ(d,ch,1)+Offset(1),Contact.XYZ(d,ch,2)+Offset(2),Contact.XYZ(d,ch,3),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Contact.NormRad(CellIndx),Fig.SphereN);
                end
                Contact.Handles(CellIndx) = surface(X,Y,Z,repmat(Contact.ColorVals(CellIndx),size(Z)),'EdgeColor','none','CDataMapping','scaled','FaceAlpha',Contact.NormAlpha(CellIndx));
                hold on;
            end
        end
    end
    set(Contact.Handles,'ButtonDownFcn',@ContactSelect);
    Fig.cbh.Label.String = Contact.DataVariable; 
end

%==================== UPDATE PLOT BASED ON GUI INPUT ======================
function ContactInput(hObj, Event, Indx)
global Contact Fig
    switch Indx
        case 1      %=========== SELECTED SESSION DATE
            Contact.Current.DateIndx = get(hObj,'Value');                                   	% Get new date selection
            Contact.Indx = min(find(Contact.CellIndxData(:,2) == Contact.Current.DateIndx));
            Contact.Current.ChannelIndx = Contact.CellIndxData(Contact.Indx,3);
            Contact.Current.CellIndx = Contact.CellIndxData(Contact.Indx,4);
            UpdateUIFields(Contact.Indx);
            UpdateHighlight(Contact);

        case 5      %=========== SELECT CHANNEL NUMBER
            ChannelIndx = str2double(get(hObj,'String'));                                      % Get new contact number
            if ~ismember(ChannelIndx,1:size(Contact.XYZ, 2))
                set(hObj,'String', Contact.Current.ChannelIndx);
            else
                Contact.Current.ChannelIndx = ChannelIndx;
                Contact.Current.DateIndx = get(Fig.Handles.CellInput(1),'value');
                AllDates = find(Contact.CellIndxData(:,2) == Contact.Current.DateIndx);             % Get indices for currently selected date
                Contact.Indx = min(AllDates(find(Contact.CellIndxData(AllDates,3)==Contact.Current.ChannelIndx))); % Find cell index
                UpdateUIFields(Contact.Indx);
                UpdateHighlight(Contact);
            end
                
        case 6      %=========== SELECT CELL NUMBER
            CellIndx = str2double(get(hObj,'String'));
            [temp, AvailableCellIndx] = ismember([Contact.Current.DateIndx, Contact.Current.ChannelIndx],Contact.CellIndxData(:,2:3), 'rows');
            if ~ismember(CellIndx, 1:max(Contact.CellIndxData(AvailableCellIndx,4)))
                set(hObj,'String', Contact.Current.CellIndx);
            else
                Contact.Current.CellIndx = CellIndx;
                Contact.Indx = find(Contact.CellIndxData(:,3)==[Contact.Current.DateIndx,Contact.Current.ChannelIndx,Contact.Current.CellIndx]);
                UpdateUIFields(Contact.Indx);
                UpdateHighlight(Contact);
            end
            
        case 7      %=========== CHANGE COLOR VALUE FOR CURRENT DATA POINT
            Contact.ColorVals(Contact.Indx) = str2double(get(hObj,'String'));                   % Get new color value
            ZData = get(Contact.Handles(Contact.Indx),'ZData');                                 % Get z-data for all vertices of current object
            NewColorData = repmat(Contact.ColorVals(Contact.Indx), size(ZData));                % Create c-data for all vertices
            set(Contact.Handles(Contact.Indx),'CData', NewColorData);                           % Apply value to render
            
         case 8      %=========== CHANGE ALPHA VALUE FOR CURRENT DATA POINT
            Contact.Alpha(Contact.Indx) = str2double(get(hObj,'String'));                     	% Get new color value
            set(Contact.Handles(Contact.Indx), 'FaceAlpha', Contact.Alpha(Contact.Indx));       % Apply value to render
                
        case 9      %=========== CHANGE RADIUS VALUE FOR CURRENT DATA POINT
            p = Contact.Indx;
        	Contact.rad(Contact.Indx) = str2double(get(hObj,'String'));                       	% Get new radius value
            RadRange = [min(Contact.rad), max(Contact.rad)];                                    % Get new range of radius values
            if numel(unique(Contact.rad))>1
                Contact.NormRad = (Contact.rad-RadRange(1))/ diff(RadRange);                  	% Normalize radius values (0-1)
                Contact.NormRad = Fig.SphereRadii(1)+(Contact.NormRad*diff(Fig.SphereRadii));  	% Scale radius values to specified range
            else
                Contact.NormRad = ones(size(Contact.rad))*Fig.SphereRadii(2);                   % If all radii supplied are same, use max
            end
            
            delete(Contact.Handles(Contact.Indx));                                           	% delete current data point
            Coordinates = Contact.XYZ(Contact.Current.DateIndx, Contact.Current.ChannelIndx, :);% Get coordinates of selected contact
            [X,Y,Z] = ellipsoid(Coordinates(1),Coordinates(2),Coordinates(3),Contact.NormRad(p),Contact.NormRad(p),Contact.NormRad(p),Fig.SphereN);
            Contact.Handles(Contact.Indx) = surface(X,Y,Z,repmat(Contact.ColorVals(Contact.Indx),size(Z)),...
                'EdgeColor','r','FaceAlpha',Contact.Alpha(Contact.Indx),'ButtonDownFcn',@ContactSelect); % Apply value to render     
    end
end

%===================== UPDATE ATLAS STRUCTURES
function StructView(hObj, Evnt, Indx)
    global Fig Structures
    switch Indx
        case 1                          %================== Change selected stucture
            Structures.CurrentStructure = get(hObj,'Value');
            set(Fig.Handles.StructInput(2), 'BackgroundColor', Structures.Colors(Structures.CurrentStructure,:));
            set(Fig.Handles.StructInput(3), 'value', Structures.Opacity(Structures.CurrentStructure)); 
            set(Fig.Handles.StructInput(4), 'value', Structures.Smoothing(Structures.CurrentStructure));
            set(Fig.Handles.StructLabel(3),'String',sprintf('Opacity = %d %%', Structures.Opacity(Structures.CurrentStructure)*100));
            set(Fig.Handles.StructInput(5),'value',Structures.On(Structures.CurrentStructure));
            set(Fig.Handles.StructInput(6),'value',Structures.Wire(Structures.CurrentStructure));
      
        case 2                          %================== Change current stuctures color
            Structures.Colors(Structures.CurrentStructure,:) = get(hObj,'BackgroundColor');
            Structures.Colors(Structures.CurrentStructure,:) = uisetcolor(Structures.Colors(Structures.CurrentStructure,:));
            set(Fig.Handles.StructInput(2), 'BackgroundColor', Structures.Colors(Structures.CurrentStructure,:));
            if Structures.Wire(Structures.CurrentStructure) == 1
                set(Structures.Handles{Structures.CurrentStructure},'edgecolor',Structures.Colors(Structures.CurrentStructure,:));
            elseif Structures.Wire(Structures.CurrentStructure) == 0
            	set(Structures.Handles{Structures.CurrentStructure},'facecolor',Structures.Colors(Structures.CurrentStructure,:));
            end
            
        case 3                          %============== Update layer transparency
            Structures.Opacity(Structures.CurrentStructure) = get(hObj,'Value');
            ValueString = sprintf('Opacity = %.0f %%', Structures.Opacity(Structures.CurrentStructure)*100);
            set(Fig.Handles.StructLabel(3),'String',ValueString);
            set(Structures.Handles{Structures.CurrentStructure},'facealpha',Structures.Opacity(Structures.CurrentStructure),'edgealpha',Structures.Opacity(Structures.CurrentStructure));

        case 4                    	%============== Update [continuous variable] ???
            Structures.Smoothing(Structures.CurrentStructure) = get(hObj,'Value');
            ValueString = sprintf('Smoothing = %.0f %%', Structures.Smoothing(Structures.CurrentStructure)*100);
            set(Fig.Handles.StructLabel(4),'String',ValueString);
%             set(Structures.Handles{Structures.CurrentStructure},
            
        case 5                      %============== Toggle structure visibility on/off
            Structures.On(Structures.CurrentStructure) = get(hObj,'Value');
            if Structures.On(Structures.CurrentStructure) == 1                             	% If layer was turned on...
                set(Structures.Handles{Structures.CurrentStructure},'visible','on');      	% make visible
            elseif Structures.On(Structures.CurrentStructure) == 0                       	% If layer was turned off...
                set(Structures.Handles{Structures.CurrentStructure},'visible','off');       % make invisible
            end
            
        case 6                      %============== Toggle mesh faces on/off
            Structures.Wire(Structures.CurrentStructure) = get(hObj,'Value');
            if Structures.Wire(Structures.CurrentStructure) == 1                        	
                set(Structures.Handles{Structures.CurrentStructure},'facecolor','none');     
                set(Structures.Handles{Structures.CurrentStructure},'edgecolor',Structures.Colors(Structures.CurrentStructure,:));
                
            elseif Structures.Wire(Structures.CurrentStructure) == 0                    
                set(Structures.Handles{Structures.CurrentStructure},'facecolor',Structures.Colors(Structures.CurrentStructure,:));      	
                set(Structures.Handles{Structures.CurrentStructure},'edgecolor','none');
            end
    end

end

%========================= DATA VIEW PANNEL ===============================
function DataView(hObj, Event, Indx)
global Contact Fig Structures Data
    switch Indx
        case 1      %==================== LOAD NEW DATA SET

            %===== DELETE CURRENT DATA
          	if isfield(Contact,'Handles')
                delete(Contact.Handles);                                                      	% Delete existing data points
            end
            delete(Fig.Data.Spatial.ph);                                                        % Delete previous plot data
            delete(Fig.Data.Hist);                                                              % Delete previous histogram data
          	if isfield(Fig,'Highlight') && ishandle(Fig.Highlight)                              % If highlight markers were plotted...
                delete(Fig.Highlight);                                                          % Delete
                delete(Fig.Data.DistHandle);
            end
%             Fig.Data.Spatial = rmfield(Fig.Data.Spatial, 'ph');                                 % Delete data handles
            clearvars -global Contact;                                                       	% Clear Contact structure from previous session
            
            %===== LOAD NEW DATA
            Data.Selected = get(hObj,'Value');
            h = msgbox(sprintf('Loading data set: %s', Data.Files{Data.Selected}),'Loading data...');
            load(Data.Files{Data.Selected});
            RadRange = [min(Contact.rad), max(Contact.rad)];                                    % Get range of radius values
            if numel(unique(Contact.rad))>1 
                Contact.NormRad = (Contact.rad-RadRange(1))/ diff(RadRange);                  	% Normalize radius values (0-1)
                Contact.NormRad = Fig.SphereRadii(1)+(Contact.NormRad*diff(Fig.SphereRadii));  	% Scale radius values to specified range
            else
                Contact.NormRad = ones(size(Contact.rad))*Fig.SphereRadii(2);                   % If all radii supplied are same, use max
            end
            AlphaRange = [min(Contact.Alpha), max(Contact.Alpha)];                              % Get range of alpha values
            % Contact.NormAlpha = (Contact.Alpha-AlphaRange(1))/ diff(AlphaRange);                % Normalize radius values (0-1)
            if diff(AlphaRange)==0
                Contact.NormAlpha = ones(size(Contact.Alpha));
            end
            LoadNewData(Contact);
            close(h);
          
        case 4      %==================== CHANGE DATA PLOT AXES
        	Fig.AxisSelected = get(hObj,'Value');
            axes(Fig.Handles.DataAx(2));    
            Allpos = reshape(Contact.XYZ(:,:,Fig.AxisSelected)',[size(Contact.XYZ,1)*size(Contact.XYZ,2), 1]);
            for i = 1:size(Contact.CellIndxData,1)
                PosIndx(i) = ((Contact.CellIndxData(i,2)-1)*size(Contact.XYZ,2))+Contact.CellIndxData(i,3);
                set(Fig.Data.Spatial.ph(i), 'xdata', Allpos(PosIndx(i)));
            end
            xlabel([Fig.Data.AxisLabels{Fig.AxisSelected},' (mm)']);
            AllLims = [Fig.Xlim; Fig.Ylim; Fig.Zlim];
            set(Fig.Handles.DataAx(2),'xlim', AllLims(Fig.AxisSelected,:));
            
        case {5, 8}	%==================== CHANGE DATA THRESHOLDS
            Thresh = [str2double(get(Fig.Handles.DataInput(5),'string')), str2double(get(Fig.Handles.DataInput(8),'string'))];
            set(Fig.Handles.DataAx(1), 'xlim', Thresh);
            set(Fig.Handles.DataAx(2), 'ylim', Thresh);
            axes(Fig.Handles.MainAx);
            NewCdata = Contact.ColorVals;
            NewCdata(NewCdata<Thresh(1)) = Thresh(1);
            NewCdata(NewCdata>Thresh(2)) = Thresh(2);
            for i = 1:numel(Contact.Handles)                    % For each data point...
                set(Contact.Handles(i),'cdata', repmat(NewCdata(i), size(get(Contact.Handles(1),'Zdata'))));
            end
            if isfield(Fig,'Highlight') && ishandle(Fig.Highlight)
                set(Fig.Highlight, 'cdata', repmat(mean(NewCdata), size(get(Fig.Highlight, 'ZData'))));
            end
            
    	case {6, 9}	%==================== CHANGE ALPHA THRESHOLDS
            AlphaThresh = [str2double(get(Fig.Handles.DataInput(6),'string')), str2double(get(Fig.Handles.DataInput(9),'string'))];
            TempAlpha = Contact.Alpha;                                                          % Copy raw alpha values to temp variable
            TempAlpha(TempAlpha < AlphaThresh(1)) = AlphaThresh(1);                             % Threshold the alpha values
            TempAlpha(TempAlpha > AlphaThresh(2)) = AlphaThresh(2);
            Contact.NormAlpha = (TempAlpha-AlphaThresh(1))/ diff(AlphaThresh);                  % Normalize alpha values (0-1)                                     
            for i = 1:numel(Contact.Handles)                                                    % For each data point...
                set(Contact.Handles(i),'facealpha', Contact.NormAlpha(i));                      % Update scaled alpha
            end
            
        case 7  	%==================== INVERT ALPHA VALUES
            Fig.Data.InvertAlpha = get(hObj,'Value');
            Contact.Alpha = 1-Contact.Alpha; 
            DataView([], [], 6);                                                                % Threshold alpha
    end

end

%====================== MAIN OPTIONS PANNEL ===============================
function OptionsInput(hObj, Event, Indx)
global Contact Fig Structures
    switch Indx
        case 1      %====================== MANUAL ROTATE
            if get(hObj,'Value')==1
                rotate3d;
            elseif get(hObj,'Value')==0
                rotate3d off;
            end
        case 2      %====================== SET CAMERA POSITION
            axes(Fig.Handles.MainAx);
            [az,el] = view;     
            ans = inputdlg({'Azimuth (deg)','Elevation (deg)'},'Set camera view',1,{num2str(az),num2str(el)});
            if ~isempty(ans)
                view(str2double(ans{1}), str2double(ans{2}));
            end
            
        case 3      %====================== SET LIGHTING POSITION TO HEADLIGHT
            axes(Fig.Handles.MainAx);
            if ishandle(Fig.Handles.Light(2))
                delete(Fig.Handles.Light(2));
            end
            Fig.Handles.Light(2)= camlight('headlight');
            
        case 4      %====================== SET AXES LIMITS
            axes(Fig.Handles.MainAx);
            Xlim = get(gca,'xlim');
            Ylim = get(gca,'ylim');
            Zlim = get(gca,'zlim');
            ans = inputdlg({'Medial-lateral (mm)','Anterior-posterior (mm)','Inferior-Superior (mm)'},'Set axes limits',1,{num2str(Xlim),num2str(Ylim),num2str(Zlim)});
            if ~isempty(ans)
                Xlim = str2double(ans{1});
                Ylim = str2double(ans{2});
                Zlim = str2double(ans{3});
                set(gca,'xlim',Xlim','ylim',Ylim,'zlim',Zlim);
            end
            
        case 5      %====================== SAVE APPEARANCE SETTINGS
            [file, path] = uiputfile('*.mat','Save appearance settings','AppearanceSettings_1.mat');
            if ischar(file)
                NewFig = Fig;                                                                       % Copy structture to new variable
                NewFig = rmfield(NewFig, 'Handles');                                              	% Remove field containing handles  
                MatFilename = fullfile(path, file);
                save(MatFilename, 'NewFig', 'Structures');
                h = msgbox(sprintf('Appearance settings saved to: %s', file),'Save sucessful!');
            end
            
        case 6      %====================== LOAD APPEARANCE SETTINGS
            [file, path] = uigetfile('*.mat','Select appearance settings');
            if ischar(file)
                MatFilename = fullfile(path, file);
                NewData = load(MatFilename);
                UpdateFigureAppearence(NewData);
            end
            
        case 7      %====================== TOGGLE ZOOM
            axes(Fig.Handles.MainAx);
            if get(hObj,'Value')==1
                Fig.ZoomH.Enable = 'on';
            elseif get(hObj,'Value')==0
                Fig.ZoomH.Enable = 'off';
            end
            
        case 8      %====================== SAVE IMAGE
            DefaultFilename = sprintf('%s_%s_N=%d.fig',Contact.Subject, Contact.DataVariable, size(Contact.ColorVals,1));
            [file, path] = uiputfile({'*.fig;*.png'}, 'Save figure as:', DefaultFilename);
            if ischar(file)
                FullFilename = fullfile(path, file);
                ImFormat =  file(findstr(file,'.')+1:end);
                set(gcf,'InvertHardcopy','off');                                
                TempFigH = figure('Color',Fig.Background,...                 	% Set the figure window background color
                                'Name','Exporting figure...',...                % Give the figure window a name
                                'Renderer','OpenGL',...                         % Use OpenGL renderer
                                'OuterPosition', Fig.Rect);                     % position figure window to fit fullscreen;
                h = copyobj(Fig.Handles.MainAx, TempFigH);                                 % Copy main 3D axes to new figure window
                colormap(Fig.Colormap);                                         % Restore current colormap
                set(h,'position',[0.1 0.1 0.8 0.8]);                            % Set axes position to fill figure
                title(sprintf('%s (N = %d)', Contact.DataVariable, size(Contact.ColorVals,1)),'horizontalalignment','left','fontsize',18);
                switch ImFormat
                    case 'png'
                        export_fig(FullFilename,['-',ImFormat]);                % Export .png image
                    case 'fig'
                        saveas(h, FullFilename);                              	% Export Matlab .fig
                end
                close(TempFigH);                                             	% Close temp figure window
            end
            
        case 9      %====================== LAUNCH OPEN GL BASED 3D VIEW
            Params = InitializeAtlasViewer3D;
            AtlasViewer3D(Params, Contacts, Structures);
    end
end


%====================== UPDATE FIGURE APPEARANCE BASED ON SAVED SETTINGS
function UpdateFigureAppearence(LoadedData)
global Fig Structure
    AllFields = fieldnames(LoadedData.NewFig)


end