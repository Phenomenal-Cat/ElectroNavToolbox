function FigHandle = GridHistory(SessionHistory, Mode, GridID)

%============================ GridHistory.m ===============================
% 
%
% INPUTS:
%   SessionHistory:  spreadsheet file (.xls/.csv) containing recording session
%               data with the following colums:
%               1) Date 
%               2) Medial-lateral grid hole coordinate
%               3) Posterior-anterior grid hole coordinate
%               4) Target depth
%               5) Guide tube length
%               6) MRI collected?
%               7) Electrode name
%   Mode:   1 = Plots number of recording sessions for each grid hole
%           2 = Plots time since last recording session for each grid hole
%           3 = Plots mean penetration depth for each grid hole
%           4 = Highlights which grid holes have MRI of electrode position
%   GridID:     string corresponding to one of the grid types defined in 
%               GetGridParams.m.
%
%==========================================================================

if nargin == 0
    Mode = 1;
    SessionHistory = 'ElectrodeLocations.xls';
end

%=========================== Set grid parameters
if ~exist('GridID','var')
    GridIDs = GetGridParams;                         	% Get list available types of grid
    [Selection,ok] = listdlg('ListString',GridIDs,'PromptString','Select grid type');
    GridID = GridIDs{Selection};
end
Grid = GetGridParams(GridID);                          	% Get grid parameters based on grid ID                 
HoleCoordinates = nan(Grid.TotalHoles,2);                   
HoleNumber = 1;
for i = 1:Grid.HolesPerDim
    for h = 1:Grid.HolesPerColumn(i)
      	x = ((-((Grid.HolesPerColumn(i)-1)/2)+(h-1))*Grid.InterHoleSpacing);
      	y = ((((Grid.HolesPerDim+1)/2)-i)*Grid.InterHoleSpacing);
        HoleCoordinates(HoleNumber,:) = [x,y];
        HoleNumber = HoleNumber+1;
    end
end


%========================== Load recording history data
[a,b,Ext] = fileparts(SessionHistory);
switch Ext
    case '.xls'
     	DateSheet = 1;
        [num,txt,raw] =  xlsread(SessionHistory,DateSheet,'');        	% Read data from Excel file
        Headers = txt{1,:};                                             % Skip row containing column titles
        num(1,:) = [];                                                  % Remove nans
        Dates = num(:,1)+datenum('30-Dec-1899');                        % Convert Excel dates to Matlab dates
    case '.csv'
        fid = fopen(SessionHistory,'rt');
        Headers = textscan(fid, '%s%s%s%s%s%s%s\n', 1, 'delimiter', ',');
        data = textscan(fid, '%f %f %f %f %f %f %s','headerlines',1,'delimiter',',');
        fclose(fid);
        Dates = data{1};   
        num(:,2:3) = cell2mat(data(:,2:3));
        num(:,4) = cell2mat(data(:,4));
end
DateStrings = datestr(Dates);  
[Selection,ok] = listdlg('ListString',DateStrings,'SelectionMode','multi','PromptString','Select sessions to include:');
if ok==0
    FigHandle = [];
    return;
end
Target = num(Selection,[2,3]);
% Target(:,2) = -Target(:,2);                                     % Invert AP axis for 2013 recording sessions
TargetDepth = num(Selection,4);                                 
Session.Date = datestr(Dates(Selection));  
DaysElapsed = datenum(date)-Dates(Selection);                   % Days since each recording


%========================== Set visualization parameters
if Mode == 1
    ColorTitle = 'Number of sessions';
    Grid.HoleTally = zeros(1, Grid.TotalHoles);
    for h = 1:numel(Target(:,1))
        HoleNumber = find(ismember(HoleCoordinates, Target(h,:), 'rows'));  % Get grid hole index
        Grid.HoleTally(HoleNumber) = Grid.HoleTally(HoleNumber)+1;          % Increase session tally count
    end
    YTicks = 0:1:max(Grid.HoleTally);
    YLims = [0 max(Grid.HoleTally)];
    ColormapMode = [Grid.RGB*0.7; cool];
    
elseif Mode == 2
    ColorTitle = 'Days since last recording';
    Grid.HoleTally = inf(1, Grid.TotalHoles);
    for h = 1:numel(Target(:,1))
        HoleNumber = find(ismember(HoleCoordinates, Target(h,:), 'rows')); 
        
        if DaysElapsed(h) < Grid.HoleTally(HoleNumber)
            Grid.HoleTally(HoleNumber) = DaysElapsed(h);
        end
    end
    Grid.HoleTally(Grid.HoleTally==inf) = 30;
    YTicks = 0:7:(7*4);
    YLims = [0 7*4];
    ColormapMode = [cool; Grid.RGB*0.7];
end


%========================= Plot data
FigBackground = [0.7 0.7 0.7];
FigHandle = figure( 'Name','Recording History',...          % Open a figure window with specified title
                    'NumberTitle','off',...                 % Remove figure number from title
                    'Color',FigBackground,...               % Set the figure window background color
                    'Menu','none','Toolbar','none');       	% Turn off toolbars to save space
GridObject = DrawGrid(Grid);
title(sprintf('Total recordings: %s - %s',DateStrings(Selection(1),:),DateStrings(Selection(end),:)),'FontWeight','bold','FontSize',14);
set(gca,'XTick',-8:2:8);
set(gca,'YTick',-8:2:8);
set(gca,'fontsize',12);
set(gca,'color',FigBackground);
Labels(1) = xlabel('Medial-Lateral');                                        
Labels(2) = ylabel('Posterior-Anterior');
grid on;
axis equal tight;
cbr = colorbar;
colormap(ColormapMode);
set(cbr,'YTick',YTicks,'YLim',YLims);
set(get(cbr,'ylabel'),'String', ColorTitle,'FontSize',12);
end


%% ========================= SUBFUNCTIONS =================================
function h = PlotCircle(x,y,r,c)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,c);
end

function h = FillCircle(target,r,N,c)
THETA=linspace(0,2*pi,N);
RHO=ones(1,N)*r;
[X,Y] = pol2cart(THETA,RHO);
X=X+target(1);
Y=Y+target(2);
h=fill(X,Y,c,'EdgeColor','none');
end

%============================== DRAW GRID =================================
function GridObject = DrawGrid(Grid)
    GridObject(1) = FillCircle([0 0],Grid.OuterRadius,100,'y');
    hold on;    
    holeno = 1;
    for i = 1:Grid.HolesPerDim
        for h = 1:Grid.HolesPerColumn(i)
            x = (-((Grid.HolesPerColumn(i)-1)/2)+(h-1))*Grid.InterHoleSpacing;
            y = (((Grid.HolesPerDim+1)/2)-i)*Grid.InterHoleSpacing;
            r = Grid.HoleDiameter/2;
            n = 20;
            c = repmat(Grid.HoleTally(holeno),[1,n]);
            GridObject(end+1) = FillCircle([x,y],r,n,c);
            holeno = holeno+1;
        end
    end
    GridObject(end+1) = plot([0 0],[-Grid.OuterRadius,Grid.OuterRadius],'-k');
    GridObject(end+1) = plot([-Grid.OuterRadius,Grid.OuterRadius],[0 0],'-k');
    GridParent = hgtransform('Parent',gca);
    set(GridObject,'Parent',GridParent);
end