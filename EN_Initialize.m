%============================ EN_Initialize.m ========================
% This function allows the user to set the default directory and file paths
% for files required by ElectroNav. These parameters are saved to a .mat
% file in the ElectroNav directory and loaded each time the GUI is run.
%
% REVISIONS:
%   03/04/2014 - Written by APM
%   30/01/2015 - Updated
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy
% � Copyleft 2014, GNU General Public License
%==========================================================================

function ParamsOut = EN_Initialize(DefaultParamsFile, SubjectID)

persistent Params 
persistent Fig
addpath(genpath('ENSubfunctions'));

%======================= CHECK INPUT
if nargin ==0 || ~exist(DefaultParamsFile,'file')
    if nargin ==0
        Msg = sprintf('A parameters file must be provided!');
    elseif ~exist(DefaultParamsFile,'file')
        Msg = sprintf('The parameters file ''%s'' was not found!',DefaultParamsFile);
    end
    Choice = questdlg(Msg,'File not found!','Select file','Continue','Cancel','Cancel');
    switch Choice
        case 'Select file'
            [filename, pathname, filterindex] = uigetfile({'*.mat','Default parameters (*.mat)'},...
                                                          'Select default parameters file');
            if filename==0
                return;
            end
            DefaultParamsFile = fullfile(pathname, filename);
        case 'Cancel'
            return;
        otherwise
            return;
    end
end

%======================= LOAD PARAMATERS FILE
load(DefaultParamsFile);
Subjects{1} = '';
for s = 1:numel(Defaults)
    Subjects{s+1} = Defaults(s).SubjectID;
end
if exist('SubjectID','var')
    SubjectIndx = find(~cellfun(@isempty, strfind(Subjects,SubjectID)))-1;
    if ~isempty(SubjectIndx)
        Params.SubjectID = Subjects{SubjectIndx+1};
    else
        error('Subject ''%s'' was not found in paramaters file ''%s''!', SubjectID, DefaultParamsFile);
    end
end
Subjects{end+1} = '*Add new subject*';                                	% Provide option to add new subject
Params.Defaults = Defaults;                                             
Params.RootDir = fileparts(mfilename('fullpath'));                      % Get ElectroNav root path




%========================= OPEN GUI WINDOW ================================
Fig.Handle = figure;%(typecast(uint8('ENav'),'uint32'));                  % Assign GUI arbitrary integer        
if strcmp('ElectroNavInit', get(Fig.Handle, 'Tag')), return; end        % If figure already exists, return
Fig.FontSize = 12;
Fig.TitleFontSize = 14;
Fig.PanelYdim = 130;
Fig.Rect = [400 400 600 600];                                           % Specify figure window rectangle
set(Fig.Handle,     'Name','ElectroNav: Initialization',...             % Open a figure window with specified title
                    'Tag','ElectroNavInit',...                          % Set figure tag
                    'Renderer','OpenGL',...                             % Use OpenGL renderer
                    'OuterPosition', Fig.Rect,...                       % position figure window
                    'NumberTitle','off',...                             % Remove figure number from title
                    'Resize', 'off',...                                 % Prevent resizing of GUI window
                    'Menu','none',...                                   % Turn off memu
                    'Toolbar','none');                                  % Turn off toolbars to save space
Fig.Background = get(Fig.Handle, 'Color');                              % Get default figure background color
Fig.FontSize = 10;
Fig.Margin = 20;   
Fig.AcceptedColor = [0 1 0];
Fig.RejectedColor = [1 0 0];
Fig.InputTags = {'HistoryFile','ExpDir','GridID','MRIanat',...        	% Declare tags for GUI input fields
                 'MRIatlas','MRIxform','MRIVTKdir','MRIElecdir'};
Fig.Fields = {'HistoryFile','ExpDir','GridID','MRI','Atlas','Xform','VTKdir','ElecMRI'};    % Set parameter field names
BoxPos{1} = [Fig.Margin,330,Fig.Rect(3)-Fig.Margin*2, 160];           	% Set group controls positions
BoxPos{2} = [Fig.Margin,70,Fig.Rect(3)-Fig.Margin*2, 240];


Logo= imread(fullfile('Documentation','ElectroNavLogo1.png'),'BackgroundColor',Fig.Background);
LogoAx = axes('box','off','units','pixels','position', [320, 520, 260, 42],'color',Fig.Background);
image(Logo);
axis off

%=================== SUBJECT SELECTION BOX
TipStr = 'Select an existing subject';
uicontrol(  'Style', 'text',...
            'String','Select subject:',...
            'Background',Fig.Background, ...
            'units','pixels',...
            'Position', [Fig.Margin*2,540,150,20],...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left');
Fig.NameH = uicontrol(  'Style', 'popup',...
            'Background', 'white',...
            'Tag', 'Subject', ...
            'String', Subjects,...
            'value', 1,...
            'units','pixels',...
            'Position', [140,540,150,20],...
            'TooltipString', TipStr,...
            'Callback', @GetSubjectData);


%================= PHYSIOLOGY DATA PANEL
Fig.Physio.Handle = uipanel( 'Title','Electrophysiology Data',...
                'FontSize',Fig.FontSize,...
                'BackgroundColor',Fig.Background,...
                'Units','pixels',...
                'Position',BoxPos{1},...
                'Parent',Fig.Handle);        
        
%================== SESSION HISTORY SPREADSHEET SELECTION
TipStr = 'Select the spreadsheet (*.xls/*.csv) containing records of previous physiology recording sessions';
uicontrol(  'Style', 'text',...
            'String','Select session history:',...
            'Background',Fig.Background, ...
            'units','pixels',...
            'Position', [Fig.Margin,100,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'Parent',Fig.Physio.Handle,...
            'HorizontalAlignment', 'left');
Fig.InH(1) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'HistoryFile', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,100,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.Physio.Handle,...
            'Callback', {@FileCheck, 1});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,100,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseHistory', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.Physio.Handle,...
            'Callback', {@FileSelect, 1});
  

%================= ANALYSED PHYSIOLOGY DATA DIRECTORY
TipStr = 'Select the directory containing analysed experimental data';
uicontrol(  'Style', 'text',...
            'String','Select experiment directory:',...
            'Background',Fig.Background, ...
            'units','pixels',...
            'Position', [Fig.Margin,60,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'Parent',Fig.Physio.Handle,...
            'HorizontalAlignment', 'left');
Fig.InH(2) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'ExpDir', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,60,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.Physio.Handle,...
            'Callback', {@FileCheck, 2});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,60,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseExp', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.Physio.Handle,...
            'Callback', {@FileSelect, 2});

        
%=================== SELECT GRID TYPE
Fig.GridIDs = GetGridParams;
Fig.GridIDs = {'',Fig.GridIDs{:}};
TipStr = 'Select the style of dural recording chamber/ grid implant';
uicontrol(  'Style', 'text',...
            'String','Select implant type:',...
            'Background',Fig.Background, ...
            'units','pixels',...
            'Position', [Fig.Margin,20,180,20],...
            'FontSize', Fig.FontSize, ...
            'Parent',Fig.Physio.Handle,...
            'HorizontalAlignment', 'left');
Fig.InH(3) = uicontrol(  'Style', 'popup',...
            'Background', 'white',...
            'Tag', 'GridID', ...
            'String', Fig.GridIDs,...
            'units','pixels',...
            'Position', [200,20,200,20],...
            'Parent',Fig.Physio.Handle,...
            'TooltipString', TipStr,...
            'Callback', @GetGridType);

        

%=========================== MRI DATA FILES ===============================
Fig.MRI.Handle = uipanel('Title','MRI Data',...
                    'FontSize',Fig.FontSize,...
                    'BackgroundColor',Fig.Background,...
                    'Units','pixels',...
                    'Position',BoxPos{2},...
                    'Parent',Fig.Handle);
                
%=================== SELECT SUBJECT'S ANATOMICAL MRI
TipStr = 'Select the anatomical (T1-weighted) MRI Nifti file (*.nii) for this subject';
uicontrol(  'Style', 'text',...
            'String','Select anatomical (.nii):',...
            'Background',Fig.Background, ...
            'parent', Fig.MRI.Handle,...
            'units','pixels',...
            'Position', [Fig.Margin,180,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left');               
Fig.InH(4) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'MRIanat', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,180,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileCheck, 4});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,180,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseAnat', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileSelect, 4});            
        
        
%=================== SELECT ALIGNED ATLAS
TipStr = 'Select the atlas Nifti file (*.nii) spatially normalized for this subject';
uicontrol(  'Style', 'text',...
            'String','Select atlas (.nii):',...
            'Background',Fig.Background, ...
            'parent', Fig.MRI.Handle,...
            'units','pixels',...
            'Position', [Fig.Margin,140,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left');               
Fig.InH(5) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'MRIatlas', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,140,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileCheck, 5});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,140,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseAtlas', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileSelect, 5});    
        
        
%=================== SELECT XFORM MATRIX FILE
TipStr = 'Select the xform transformation matrix file (*.xform) for this subject''s grid';
uicontrol(  'Style', 'text',...
            'String','Select transform (.xform):',...
            'Background',Fig.Background, ...
            'parent', Fig.MRI.Handle,...
            'units','pixels',...
            'Position', [Fig.Margin,100,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left');               
Fig.InH(6) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'MRIxform', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,100,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileCheck, 6});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,100,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseXform', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileSelect, 6});    

        
%=================== SELECT VTK DIRECTORY
TipStr = 'Select the directory containing 3D mesh files (*.vtk) for this subject''s anatomical structures';
uicontrol(  'Style', 'text',...
            'String','Select 3D mesh directory:',...
            'Background',Fig.Background, ...
            'parent', Fig.MRI.Handle,...
            'units','pixels',...
            'Position', [Fig.Margin,60,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left');               
Fig.InH(7) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'MRIVTKdir', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,60,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileCheck, 7});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,60,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseVTK', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileSelect, 7});      
        
        %=================== SELECT VTK DIRECTORY
TipStr = 'Select the directory containing post-recording MRI scans (*.nii) for this subject.';
uicontrol(  'Style', 'text',...
            'String','Select electrode MRI directory:',...
            'Background',Fig.Background, ...
            'parent', Fig.MRI.Handle,...
            'units','pixels',...
            'Position', [Fig.Margin,20,180,20],...
            'TooltipString', TipStr,...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left');               
Fig.InH(8) = uicontrol(  'Style', 'edit',...
            'Background', [0.3 0.3 0.3],...
            'Tag', 'MRIElecdir', ...
            'String', '',...
            'units','pixels',...
            'Position', [200,20,300,25],...
            'TooltipString', TipStr,...
            'Enable','off',...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileCheck, 8});
uicontrol(  'Style', 'pushbutton', ...
            'units','pixels',...
            'Position', [520,20,20,20], ...
            'FontSize', Fig.FontSize, ...
            'String', '...', 'Tag', 'BrowseElecMRI', ...
            'TooltipString', 'Select file', ...
            'Parent',Fig.MRI.Handle,...
            'Callback', {@FileSelect, 8}); 
        

%================ SAVE DEFAULTS FOR THIS SUBJECT?
uicontrol(  'Style', 'pushbutton',...
            'String','Continue',...
            'parent', Fig.Handle,...
            'tag','Cont',...
            'units','pixels',...
            'Position', [Fig.Margin,20,120,30],...
            'TooltipString', 'Use current inputs',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 1});   
uicontrol(  'Style', 'pushbutton',...
            'String','Save',...
            'parent', Fig.Handle,...
            'tag','Save',...
            'units','pixels',...
            'Position', [160,20,120,30],...
            'TooltipString', 'Save current inputs to file',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 2});    
uicontrol(  'Style', 'pushbutton',...
            'String','Cancel',...
            'parent', Fig.Handle,...
            'tag','Cancel',...
            'units','pixels',...
            'Position', [300,20,120,30],...
            'TooltipString', 'Exit',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 3});         
uicontrol(  'Style', 'pushbutton',...
            'String','Help',...
            'parent', Fig.Handle,...
            'tag','Help',...
            'units','pixels',...
            'Position', [440,20,120,30],...
            'TooltipString', 'Open instructions document',...
            'FontSize', Fig.FontSize, ...
            'HorizontalAlignment', 'left',...
            'Callback', {@OptionSelect, 4}); 

hs = guihandles(Fig.Handle);                                % get UI handles
guidata(Fig.Handle, hs);                                    % store handles
set(Fig.Handle, 'HandleVisibility', 'callback');            % protect from command line



%================= AUTO-FILL GUI FIELDS
if exist('SubjectIndx','var')==1
    set(Fig.NameH, 'value', SubjectIndx+1);
    FillField(SubjectIndx);
end

drawnow;
uiwait(Fig.Handle);
ParamsOut = Params;


%% ====================== GUI CALLBACK FUNCTIONS ==========================

    %==================== FILL GUI FIELDS
    function FillField(SubjectIndx)
        for f = 1:numel(Fig.Fields)                                                                             % For each field...
            if isfield(Params.Defaults, Fig.Fields{f})                                                          % If a default exists...
                eval(sprintf('Params.%s = Params.Defaults(%d).%s;', Fig.Fields{f}, SubjectIndx, Fig.Fields{f}));% Copy parameter to structure
%                 Fig.InputTags{f}
%                 temp = findobj('Tag', Fig.InputTags{f})                                                        % Get GUI field handles
%                 Handle(f) = temp(1);                                                                            % Get GUI field handle
                Handle(f) = Fig.InH(f);
                File{f} = eval(sprintf('Params.Defaults(%d).%s;', SubjectIndx, Fig.Fields{f}));                 
                if strcmpi(Fig.InputTags{f}, 'GridID')
                    GridIndx = find(strcmp(Fig.GridIDs, Params.GridID));
                    if isempty(GridIndx)
                        GridIndx = 1;
                    end
                    set(Handle(f),'value', GridIndx);
                    eval('sprintf(''Params.%s = %s;'', Fig.Fields{f}, File{f});');
                else
                    set(Handle(f), 'String',File{f});                                   % Add to GUI field
                    if exist(File{f}) ~=0                                              	% If the specified file/ directory exists...
                        set(Handle(f),'Enable','on','Background',Fig.AcceptedColor);    % Set it as accepted entry in GUI 
                        eval('sprintf(''Params.%s = %s;'', Fig.Fields{f}, File{f});');  
                    elseif exist(File{f})==0                                          	% If the specified file/ directory doesn't exist...
                        set(Handle(f),'Enable','on','Background', Fig.RejectedColor);  	% Set it as rejected entry in GUI   
                        eval(sprintf('Params.%s = [];', Fig.Fields{f}));
                    end
                end
            elseif ~isfield(Params.Defaults, Fig.Fields{f})
                eval(sprintf('Params.%s = [];', Fig.Fields{f}));   
            end
        end
    end


    %==================== GET DATA FOR SELECTED SUBJECT
    function GetSubjectData(hObj, Event, Indx)
        Subjects = get(hObj,'String');
        Params.SubjectID = Subjects{get(hObj,'Value')};
        switch Params.SubjectID
            case '*Add new subject*'    %========================== Create new subject
                Subject = inputdlg('Subject ID:', 'Add new subject');                   % Get new subject name
                if isempty(Subject)
                    return;
                end
                Subjects{end} = Subject{1};                                             % Add subject to GUI menu list
                Subjects{end+1} = '*Add new subject*';
                SubjectField =  findobj('Tag', 'Subject');                          
                set(SubjectField, 'String', Subjects, 'Value', numel(Subjects)-1);      % Set GUI list to new subject
                mkdir(fullfile(Params.RootDir,'Subjects',Subject{1}));                	% Create subject directory
                Params.SubjectID = Subject{1};                                         	% Save subject ID to params structure

            case ''                     %========================== Empty subject field
                for f = 1:numel(Fig.InputTags)
                    Handle(f) = findobj('Tag', Fig.InputTags{f});
                    if strcmpi(Fig.InputTags{f}, 'GridID')
                        set(Handle(f),'Value', 1);
                    else
                        set(Handle(f),'String','','Enable','off');
                    end
                end

            otherwise                   %========================== Subject was selected
                if isfield(Params,'Defaults')                                                       % If default parameters were loaded from file...
                    for i = 1:numel(Defaults)
                        Subjects{i} = Defaults(i).SubjectID;
                    end
                    SubjectIndx = find(~cellfun(@isempty, strfind(Subjects,Params.SubjectID))); 
                    if ~isempty(SubjectIndx)                                                      	% If default parameters exist for selected subject... 
                        FillField(SubjectIndx);                                                     % Load parameters
                    end
                    Params.SubjectID = Subjects{SubjectIndx};
                end
        end

    end

    %==================== SELECT FILE/ DIRECTORY
    function FileSelect(hObj, Event, Indx)
        Fig.SelectDir = [0 1 0 0 0 0 1 1];
        Fig.FileTypes = {'*.xls;*.csv','','','*.nii;*.hdr;*.img','*.nii;*.hdr;*.img','*.xform;*.mat','',''};
        Fig.SelectPrompt = {'History file','experiment directory','','anatomical volume','atlas volume','xform matrix','VTK directory','post-session MRI directory'};
        DefaultDir = fullfile(cd, 'Subjects', Params.SubjectID);
        if Fig.SelectDir(Indx)==1
            Selection = uigetdir(DefaultDir, ['Select ',Fig.SelectPrompt{Indx}]);
            if Selection == 0
                return;
            end
        elseif Fig.SelectDir(Indx)==0
            [Filename, Pathname, FileIndx] = uigetfile(Fig.FileTypes{Indx}, ['Select ',Fig.SelectPrompt{Indx}], DefaultDir); 	% Ask user to specify file
            if Filename == 0
                return;
            end
            Selection = fullfile(Pathname, Filename);
        end
        Handle = findobj('Tag', Fig.InputTags{Indx});
        set(Handle, 'String', Selection, 'Background', Fig.AcceptedColor,'Enable', 'on');
        eval(sprintf('Params.%s = Selection;', Fig.Fields{Indx}));
    end

    %==================== CHECK FILE/ DIRECTORY INPUT
    function FileCheck(hObj, Event, Indx)
        InputString = get(hObj,'String');
        Handle = findobj('Tag', Fig.InputTags{Indx});
        if Fig.SelectDir(Indx)==1
            Type = 'dir';
        elseif Fig.SelectDir(Indx)==0
            Type = 'file';
        end
        if exist(InputString, Type)~=0
            set(Handle, 'String', InputString, 'Background', Fig.AcceptedColor,'Enable', 'on');
            eval(sprintf('Params.%s = InputString;', Fig.Fields{Indx}));
        elseif exist(InputString, Type)==0
            set(Handle, 'String', InputString, 'Background', Fig.RejectedColor);
            eval(sprintf('Params.%s = [];', Fig.Fields{Indx}));
        end
    end

    %==================== GET GRID TYPE
    function GetGridType(hObj, Event, Indx)
        GridIDs = get(hObj,'String');
        Params.GridID = GridIDs{get(hObj,'Value')};

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
                [t, CompName] = system('hostname');
                DefaultFilename = sprintf('ParamsFile_%s.mat', CompName(1:end-1));
                [Filename, Pathname, Indx] = uiputfile('*.mat','Save default parameters file',DefaultFilename);
                if Filename == 0
                    return;
                end
                save(fullfile(Pathname, Filename),'Defaults');
                msgbox('Default parameters saved!','Saved');

            case 3      %================ Cancel
                close(Fig.Handle);      % Close GUI figure
                ParamsOut = [];
                error('User cancelled ElectroNav initialization!');

        end
    end

end