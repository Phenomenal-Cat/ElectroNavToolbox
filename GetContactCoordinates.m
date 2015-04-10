function [Contact] = GetContactCoordinates(HistoryFile, Contact, Xform)

%========================= GetContactCoordinates.m ========================
%
%
% INPUTS:
%   HistoryFile:    Full path of an Excel (.xls) or Matlab data (.mat) file 
%                   containing a history of previous recording sessions. 
%   Xform:          Optional 4 x 4 transformation matrix for converting
%                   output coordinates from 'Grid space' into a native
%                   space or atlas space.
%
%   Contact.Dates:      A 1 x d cell array of date strings of recording
%                       session(s), in any standard Matlab date format.
%   Contact.DateFormat: Format off date strings provided e.g.'yyyymmdd';
%   Contact.ColorVals: 	A d x n matrix containing a scalar value to specify the 
%                       color for each data point via a colormap.
%   Contact.rad:        A d x n matrix containing a scalar value to specify the
%                       radius of each data point in millimetres.
%   Contact.Alpha:      A d x n matrix containing a scalar value to specify the
%                       opacity (1-transparency) of each data point.
%
% OUTPUTS:
%   Contact.xyz:   	1 x d array of n x 3 matrices, where d is the number of 
%                   dates input and n is the number of contacts recorded from
%                   for each session. Each matrix contains Cartesian 3D 
%                   coordinates in grid coordinate space.
%   Contact.xyz:   	1 x d array of n x 3 matrices, where d is the number of 
%                   dates input and n is the number of contacts recorded from
%                   for each session. Each matrix contains Cartesian 3D 
%                   coordinates in atlas space (transformed by Xform if supplied). 
%
% EXAMPLE:
%   HistoryFile = '/Subjects/Layla/LaylaElectrodeLocations.xls';
%   Contact.Dates = {'27-Mar-2014','15-Apr-2014','22-Apr-2014','24-Apr-2014'};
%   Contact.DateFormat = 'DD-MMM-YYYY';
%   ContactsToHightlight = {[10,11], [18,19,24],[20,21],[1,2,7]};
%   Contact.ColorVals = zeros(numel(Contact.Dates),24);
%   for d = 1:numel(Contact.Dates)
%       Contact.ColorVals(repmat(d,[1,numel(ContactsToHightlight{d})]),ContactsToHightlight{d}) = 1;
%   end
%   Contact.rad = (Contact.ColorVals+1)*0.05;
%   Contact.Alpha = (Contact.ColorVals+1)/2;
%   Contact = GetContactCoordinates(HistoryFile, Contact,[]);
%   
% REVISIONS:
%   27/04/2014 - Written by APM
%
%     ___  ______  __   __
%    /   ||  __  \|  \ |  \    APM SUBFUNCTIONS
%   / /| || |__/ /|   \|   \   Aidan P. Murphy - murphyap@mail.nih.gov
%  / __  ||  ___/ | |\   |\ \  Section of Cognitive Neurophysiology and Imaging
% /_/  |_||_|     |_| \__| \_\ Laboratory of Neuropsychology, NIMH
%==========================================================================
global Fig Contact
if nargin == 0
    HistoryFile = '/Subjects/Layla/LaylaElectrodeLocations.xls';
    SheetSelection = 1;
    Contact = struct;
end
SheetSelection = 1;

%======================== LOAD RECORDING HISTORY DATA =====================
if exist('HistoryFile','var') && exist(HistoryFile, 'file')==0                      % If the specified history file does not exist...
    if exist(fullfile(cd, HistoryFile), 'file')==2                                  % Check whether a relative path was provided...
        HistoryFile = fullfile(cd, HistoryFile);
    else
        h = warndlg(sprintf('Specified history file ''%s'' does not exist!', HistoryFile),'Save failed!','modal');      % Inform user
        uiwait(h);
        clear HistoryFile;                                                          % Delete the provided history file variable
    end
end
if exist('HistoryFile','var')==0                                                    % If history file was not specified...
    [Filename, Pathname, Indx] = uigetfile('*.xls; *.mat', 'Save current session to...'); 	% Ask user to specify Excel/ mat file to save to
    if isequal(Filename,0) || isequal(Pathname,0)                                               
        return
    end
    HistoryFile = fullfile(Pathname, Filename);                                     % Set full path of History file
end
if strcmpi(HistoryFile(end-2:end),'xls')                                            % If selected history file was Excel format...
    [status,SheetNames] = xlsfinfo(HistoryFile);                                    % Get Excel sheet names
    if ~exist('SheetSelection','var')
        [SheetSelection,ok] = listdlg('ListString',SheetNames,...                	% Ask user to select a sheet
                                 'ListSize',[160 60],...
                                 'SelectionMode', 'multiple',...
                                 'PromptString','Select Excel sheet(s):'); 
        if ok==0, return; end
    end
    Data = [];
    for s = 1:numel(SheetSelection)
        [num,txt,raw] =  xlsread(HistoryFile, SheetNames{s},'','basic');            % Read Excel file
        Data = [Data; num];
    end
    Headers = txt{1,:};
    num(1,:) = [];                                                               	% Remove nans
    DateNums = num(:,1)+datenum('30-Dec-1899');                                   	% Convert Excel dates to Matlab dates
    DateStrings = datestr(DateNums);                                
    
elseif strcmpi(Defaults.HistoryFile(end-2:end),'mat')                               % If selected file was .mat format...                                
    load(Defaults.HistoryFile);                                                     
    num = SessionHistory;
    DateStrings = datestr(num(:,1));
end

%=================== GET ELECTRODE POSITION FOR EACH DATE =================
if isfield(Contact,'DateFormat')
   Contact.DateFormat = lower(Contact.DateFormat); 
end
if isfield(Contact,'Dates')==1
    if ischar(Contact.Dates)
        Date = Contact.Dates;
        Contact = rmfield(Contact,'Dates');
        Contact.Dates{1} = Date;
    end
    for d = 1:numel(Contact.Dates)
        Contact.Dates{d} = datestr(datenum(Contact.Dates{d},Contact.DateFormat));  	% Convert dates to DD-MMM-YYYY format
        SelectionVal = strmatch(Contact.Dates{d},DateStrings);    
        if ~isempty(SelectionVal)
            Selection(d) = SelectionVal;
        else
            fprintf('WARNING: Date ''%s'' does not exist in %s!\n', Contact.Dates{d}, HistoryFile);
        end
    end              
elseif isfield(Contact,'Dates')==0                                              	% If session date/s were not specified...
    [Selection,ok] = listdlg('ListString',DateStrings,'SelectionMode','multi','PromptString','Select date of previous session:');
    if ok==0, return; end
    Contact.Dates = cellstr(datestr(DateNums(Selection)));
end
GridHole = num(Selection,[2,3]);
TipDepth = num(Selection,4);
ElectrodeID = txt(Selection+1,7);
GuideLength = num(Selection, 5);

%==================== CACLULATE POSITION OF EACH CONTACT ==================
AllElectrodeTypes = GetElectrodeParams;
xyz = [];
if exist('Contacts','var')==0 
    Contacts = cell(1,numel(Contact.Dates));
end
Contact.DatePerDatum = {};
Contact.Numbers = [];
for d = 1:numel(Contact.Dates)
    fprintf('Processing %s: electrode ID %s...\n', Contact.Dates{d}, ElectrodeID{d});
    if ~isempty(strfind(ElectrodeID{d},'_'))
        ElectrodeID{d} = ElectrodeID{d}(1:strfind(ElectrodeID{d},'_')-1);
    end
    Electrode = GetElectrodeParams(ElectrodeID{d});
    if isempty(Contacts{d})
        Contacts{d} = 1:Electrode.ContactNumber;
    end
    for c = Contacts{d}
        Contact.DatePerDatum{end+1} = Contact.Dates{d};
        Contact.Numbers(end+1) = c;
        xyz(end+1,[1,2]) = GridHole(d,:);
        xyz(end,3) = -TipDepth(d)+Electrode.TipLength+((c-1)*Electrode.ContactSpacing);
    end
end


%========================= PLOT 3D CLOUD ==================================
Fig.scnsize = get(0,'ScreenSize');                          % Get screen resolution
Fig.Rect = [0 0 Fig.scnsize(3)/2, Fig.scnsize(4)];          
Fig.FontSize = 10;
Fig.Background = repmat(0.75,[1,3]);                        
Fig.AxesBkgColor = repmat(0.75,[1,3]);
Fig.Handle = figure( 'Name','Plot Contacts',...            	% Open a figure window with specified title
                    'Color',Fig.Background,...              % Set the figure window background color
                    'Renderer','OpenGL',...               	% Use OpenGL renderer
                    'OuterPosition', Fig.Rect,...          	% position figure window to fit fullscreen
                    'NumberTitle','off',...                 % Remove figure number from title
                    'IntegerHandle','off');                 % Don't use integer handles
Fig.HighlightColor = [1 0 0];                               % Set color for highlighted contacts
Fig.HighlightAlpha = 0.3;                                   % Set alpha transparency for outer highlight volume
Fig.HighlightRad = 0.5;                                     % Set radius (mm) of outer highlight volume
MaxRadius = 0.1;                                            % Set maximum radius (mm) of contacts
N = 20;                                                     % Set number of segments per sphere

%=============== CHECK DATA FORMAT
if isfield(Contact,'ColorVals') == 0
    fprintf('No data fields supplied!');
    Contact.ColorVals = rand(1,size(xyz,1));                    
    Contact.rad = rand(1,size(xyz,1));
    Contact.Alpha = rand(1, size(xyz,1));
end
Contact.xyz = xyz;
if find(size(Contact.ColorVals)==numel(Contact.Dates))==1
    Contact.ColorVals = permute(Contact.ColorVals,[2,1]);
    Contact.rad = permute(Contact.rad,[2,1]);
    Contact.Alpha = permute(Contact.Alpha,[2,1]);
end
Contact.ColorVals = reshape(Contact.ColorVals, [1, numel(Contact.ColorVals)]);
Contact.rad = reshape(Contact.rad, [1, numel(Contact.rad)]);
Contact.Alpha = reshape(Contact.Alpha, [1, numel(Contact.Alpha)]);
if numel(find(Contact.rad>MaxRadius))>0
    fprintf('Rescaling data points to within maximum radius range!\n');
    Contact.rad = Contact.rad/max(Contact.rad(:))*MaxRadius;
end

%================ PLOT DATA
Ax(1) = axes('units','normalized','position',[0.05,0.1,0.6,0.85]);
for p = 1:size(xyz,1)
    [X,Y,Z] = ellipsoid(xyz(p,1),xyz(p,2),xyz(p,3),Contact.rad(p),Contact.rad(p),Contact.rad(p),N);
    ContactHandles(p) = surface(X,Y,Z,repmat(Contact.ColorVals(p),size(Z)),'EdgeColor','none','FaceAlpha',Contact.Alpha(p),'CDataMapping','scaled');
    hold on;
end
shading interp; 
camlight right; 
lighting phong; 
view(60,30);
grid on;
axis equal;
xlabel('Lateral-medial');
ylabel('Posterior-anterior');
zlabel('Inferior-superior');
set(gca,    'fontsize',Fig.FontSize,...
            'color',Fig.AxesBkgColor,...
            'Xlim',xlim+[-1 1],...
            'Ylim',ylim+[-1 1]);
cbh = colorbar('SouthOutside');
set(cbh, 'Position', [0.05 0.02 0.4 0.02]);
set(ContactHandles,'ButtonDownFcn',@ContactSelect);
Contact.Handles = ContactHandles;
% stereoview;


%% =========================== ADD GUI PANELS ===============================

%========================= CONTACT INFO PANEL
BoxPos = [Fig.Rect(3)-280 20 260 220];
Fig.UIhandle = uipanel( 'Title','Selected cell','FontSize',12,...
                        'BackgroundColor',Fig.Background,...
                        'Units','pixels','Position',BoxPos,...
                        'Parent',Fig.Handle);
Fig.Labels = {'Date','[x,y,z] (mm)','Contact #','Color value','Alpha value','Radius value'};
Fig.Tags = {'Date','x','y','z','ContactNo','Color', 'Alpha','Radius'};
Fig.Formats = {'%s','%.1f','%.1f','%.1f','%d','%.2f','%.2f','%.2f'};
Ypos = BoxPos(4)-20;
f = 1;
for n = 1:numel(Fig.Labels)
    Ypos = Ypos-30;
    Fig.LabelHandle(n) = uicontrol(	'Style','Text',...
                                'String',Fig.Labels{n},...
                             	'HorizontalAlignment','Left',...
                              	'pos',[10, Ypos, 100, 25],...
                                'BackgroundColor', Fig.Background,...
                              	'parent',Fig.UIhandle);
 	if n == 1
       Fig.InputHandle(f) = uicontrol(	'Style','popup',...
                            'String',Contact.Dates,...
                            'tag',Fig.Tags{f},...
                            'pos',[90, Ypos, 150, 25],...
                            'BackgroundColor', Fig.Background,...
                            'Callback',{@ContactInput,f},...
                            'parent',Fig.UIhandle);      
     	f = f+1;
 	elseif n == 2
        Xpos = 90;
        for i = 1:3
           Fig.InputHandle(f) = uicontrol(	'Style','Edit',...
                            'String','',...
                            'tag',Fig.Tags{f},...
                            'pos',[Xpos, Ypos, 48, 25],...
                            'BackgroundColor', Fig.Background,...
                            'Callback',{@ContactInput,f},...
                            'parent',Fig.UIhandle); 
            Xpos = Xpos+50;
            f = f+1;
        end
    else
        Fig.InputHandle(f) = uicontrol(	'Style','Edit',...
                                    'String','',...
                                    'tag',Fig.Tags{f},...
                                    'pos',[90, Ypos, 150, 25],...
                                    'BackgroundColor', Fig.Background,...
                                    'Callback',{@ContactInput,f},...
                                    'parent',Fig.UIhandle);
    	f = f+1;
    end
end


%======================== USER ACTIONS PANEL
BoxPos = [Fig.Rect(3)-280 260 260 220];
Fig.Actionhandle = uipanel( 'Title','Options','FontSize',12,...
                        'BackgroundColor',Fig.Background,...
                        'Units','pixels','Position',BoxPos,...
                        'Parent',Fig.Handle);
Fig.Labels = {'Manual rotate','Auto rotate','Save image','Stereo render'};
Fig.ButtonSyles = {'Togglebutton','Pushbutton','Pushbutton','Pushbutton'};
Ypos = BoxPos(4)-20;
for n = 1:numel(Fig.Labels)
    Ypos = Ypos-30;
    Fig.ButtonHandle(n) = uicontrol('Style',Fig.ButtonSyles{n},...
                                'String',Fig.Labels{n},...
                             	'HorizontalAlignment','Left',...
                              	'pos',[10, Ypos, 100, 25],...
                                'BackgroundColor', Fig.Background,...
                                'Callback',{@OptionsInput,n},...
                              	'parent',Fig.Actionhandle);
end


end


%======================== DISPLAY SELECTED CELL DATA ======================
function ContactSelect(objectHandle, eventData)
global Contact Fig
    axesHandle  = get(objectHandle,'Parent');                   % Get axes handle
    Contact.Indx = find(Contact.Handles==objectHandle);          % Get contact index number
 	Coordinates = Contact.xyz(Contact.Indx,:);                   % Get coordinates of selected contact
    
    %=========== MOVE HIGHLIGHT VOLUME IN PLOT
    set(Contact.Handles,'EdgeColor','none');                    % Turn edges off for all contacts
    set(objectHandle,'EdgeColor',Fig.HighlightColor);         	% Turn edges of selected contact red
    if isfield(Fig,'Highlight')
        if ishandle(Fig.Highlight)
            delete(Fig.Highlight);
        end
    end
    [X,Y,Z] = ellipsoid(Coordinates(1),Coordinates(2),Coordinates(3),Fig.HighlightRad,Fig.HighlightRad,Fig.HighlightRad,50);
    Fig.Highlight = surface(X,Y,Z,rand(size(Z)),'EdgeColor','none','FaceColor',Fig.HighlightColor,'FaceAlpha',Fig.HighlightAlpha);
     
    UpdateUIFields(Contact.Indx);
end

%==================== UPDATE VALUES IN UI FIELDS ==========================
function UpdateUIFields(ContactIndx)
global Contact Fig
    Contact.SelectedDateIndx = strmatch(Contact.DatePerDatum{ContactIndx},Contact.Dates);
    Coordinates = Contact.xyz(ContactIndx,:);                       % Get coordinates of selected contact
    Inputs = {Contact.DatePerDatum{ContactIndx}, Coordinates(1),Coordinates(2),Coordinates(3),...
              Contact.Numbers(ContactIndx),Contact.ColorVals(ContactIndx),Contact.Alpha(ContactIndx),Contact.rad(ContactIndx)};
    for n = 1:numel(Fig.Tags)                                       % For each GUI field tag...
        if n==1                                                     % For date popup menu...
            DateIndx = strmatch(Inputs{n},Contact.Dates)           % Get index of selected date
            set(findobj('Tag',Fig.Tags{n}), 'Value', DateIndx);
        else
            Inputs{n} = sprintf(Fig.Formats{n},Inputs{n});        	% Convert numbers to correct format string
            set(findobj('Tag',Fig.Tags{n}), 'String', Inputs{n});   % Set field to new value
        end
    end
end

%==================== UPDATE PLOT BASED ON GUI INPUT ======================
function ContactInput(hObj, Event, Indx)
global Contact Fig
    
    switch Indx
        case 1      %=========== HIGHLIGHT FIRST CONTACT FOR SELECTED DATE
            DateIndx = get(hObj,'Value');                                           % Get new date selection
            Contact.SelectedDateIndx = DateIndx;
            ContactIndices = strmatch(Contact.Dates{DateIndx},Contact.DatePerDatum);
            Contact.Indx = min(ContactIndices);
            UpdateUIFields(Contact.Indx);
            
            Coordinates = Contact.xyz(Contact.Indx,:);                              % Get coordinates of selected contact
            set(Contact.Handles,'EdgeColor','none');                                % Turn edges off for all contacts
            set(Contact.Handles(Contact.Indx),'EdgeColor',Fig.HighlightColor);      % Turn edges of selected contact red
            if isfield(Fig,'Highlight')
                if ishandle(Fig.Highlight)
                    delete(Fig.Highlight);
                end
            end
            [X,Y,Z] = ellipsoid(Coordinates(1),Coordinates(2),Coordinates(3),Fig.HighlightRad,Fig.HighlightRad,Fig.HighlightRad,50);
            Fig.Highlight = surface(X,Y,Z,rand(size(Z)),'EdgeColor','none','FaceColor',Fig.HighlightColor,'FaceAlpha',Fig.HighlightAlpha);
            
        case [2,3,4]
            
            
        case 5      %=========== MOVE HIGHLIGHT VOLUME IN PLOT
            ContactNum = str2num(get(hObj,'String'));                               % Get new contact number (0-32)
            ContactIndices = find(Contact.Numbers==ContactNum);                     % Get all possible indices
            Contact.Indx = ContactIndices(Contact.SelectedDateIndx);                % Select contact index for currently selected date
            Coordinates = Contact.xyz(Contact.Indx,:);                              % Get coordinates of selected contact
            set(Contact.Handles,'EdgeColor','none');                                % Turn edges off for all contacts
            set(Contact.Handles(Contact.Indx),'EdgeColor',Fig.HighlightColor);      % Turn edges of selected contact red
            if isfield(Fig,'Highlight')
                if ishandle(Fig.Highlight)
                    delete(Fig.Highlight);
                end
            end
            [X,Y,Z] = ellipsoid(Coordinates(1),Coordinates(2),Coordinates(3),Fig.HighlightRad,Fig.HighlightRad,Fig.HighlightRad,50);
            Fig.Highlight = surface(X,Y,Z,rand(size(Z)),'EdgeColor','none','FaceColor',Fig.HighlightColor,'FaceAlpha',Fig.HighlightAlpha);
            
        case 6      %=========== CHANGE COLOR VALUE FOR CURRENT DATA POINT
            Contact.ColorVals(Contact.Indx) = str2num(get(hObj,'String'));                  % Get new color value
            ZData = get(Contact.Handles(Contact.Indx),'ZData');                             % Get z-data for all vertices of current object
            NewColorData = repmat(Contact.ColorVals(Contact.Indx), size(ZData));            % Create c-data for all vertices
            set(Contact.Handles(Contact.Indx),'CData', NewColorData);                       % Apply value to render
            
         case 7      %=========== CHANGE ALPHA VALUE FOR CURRENT DATA POINT
            Contact.Alpha(Contact.Indx) = str2num(get(hObj,'String'));                     	% Get new color value
            set(Contact.Handles(Contact.Indx), 'FaceAlpha', Contact.Alpha(Contact.Indx)); 	% Apply value to render
                
        case 8      %=========== CHANGE RADIUS VALUE FOR CURRENT DATA POINT
            p = Contact.Indx;
        	Contact.rad(p) = str2num(get(hObj,'String'));                         	% Get new radius value
            delete(Contact.Handles(p));                                            	% delete current data point
            Coordinates = Contact.xyz(p,:);                                        	% Get coordinates of selected contact
            [X,Y,Z] = ellipsoid(Coordinates(1),Coordinates(2),Coordinates(3),Contact.rad(p),Contact.rad(p),Contact.rad(p),20);
            Contact.Handles(p) = surface(X,Y,Z,repmat(Contact.ColorVals(p),size(Z)),...
                'EdgeColor','none','FaceAlpha',Contact.Alpha(p),'ButtonDownFcn',@ContactSelect); % Apply value to render
            
    end


end


function OptionsInput(hObj, Event, Indx)
global Contact Fig
    switch Indx
        case 1      %=========== MANUAL ROTATE
            if get(hObj,'Value')==1
                rotate3d;
            elseif get(hObj,'Value')==0
                rotate3d off;
            end
        case 2      %=========== AUTO ROTATE
           	if get(hObj,'Value')==1
                rotate3d;
                [az el] = view;
                for theta = 1:5:360
                   view(az+theta, el); 
                   drawnow;
                end
                set(hObj,'Value',0);
                rotate3d off;
            end
        case 3      
            
        case 4      %=========== LAUNCH OPEN GL BASED 3D VIEW
            Params = InitializeAtlasViewer3D;
%             AtlasViewer3D(Params, Contacts, Surfaces)
    end
end