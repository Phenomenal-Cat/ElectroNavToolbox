function [Fig, Mesh] = EN_ImplantRecon(MeshDir)

%=========================== EN_ImplantRecon.m ============================
% This function reads in 3D surface mesh data (saved in .stl format) 
% from the specified directory and plots it in a figure window. Interactive
% GUI controls allow the user to manually position the recording chamber
% and record the selected trajectory in stereotaxic coordinates for use
% during surgery, .
% The directory 'MeshDir' should contain files with the following suffixes:
%
%       '*brain.stl'    segmented brain surface
%       '*ROI.stl':     manually created target region of interest surface
%       '*head.stl':    soft tissue surface reconstruction from CT scan
%       '*skull.stl':   skull surface reconstruction from CT scan
%  
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

ENpath = cd;
addpath(genpath(ENpath))

if nargin == 0
    MeshDir = uigetdir(cd, 'Select directory containing meshes to load');
end

%================ Load surface mesh data
Mesh.FileTypes   = {'head','skull','brain','ROI'};
Mesh.FileFormat  = 'stl';
for f = 1:numel(Mesh.FileTypes)
    MeshFile = wildcardsearch(MeshDir, sprintf('*%s.%s', Mesh.FileTypes{f}, Mesh.FileFormat));   % Find all .stl files in mesh directory
    if ~isempty(MeshFile)
        Mesh.Filenames{f} = MeshFile{1};
    end
end

%================ Set default appearance
SL_offset       = [0,21,15.5];              % Offset between anterior commisure (MRI-origin) and ear-bar zero (stereotaxic origin) (mm)
ChamberOffset   = [9,-8,42];                % Offset between chamber origin and ear-bar zero (stereotaxic origin) (mm)
ApplyOffset     = [1,1,1,2,2,3];            % S
theta           = -30;

ChamberTform    = [1,0,0;0,cosd(theta),sind(theta);0,-sind(theta),cosd(theta)];
TowerTform      = makehgtform('zrotate',pi/2);
TowerTform      = TowerTform(1:3,1:3);
Fig.Xlims           = [-50, 50];
Fig.Ylims           = [-40, 60];    
Fig.Zlims           = [-30, 80];

Mesh.FaceAlpha      = [0.05, 0.2, 0.2, 0.5, 0.5];                               % Set alpha/ transparency level for each mesh
Mesh.FaceColor      = [0.5 0.5 0.5; 0.8 0.8 0.8; 1 0.7 0.7; 1 0 0; 1 1 0; 1 1 0];     
Mesh.PlaneAlpha     = 0.1;
Mesh.PlaneOn        = 'on';
Mesh.ChamberOn      = 1;
Mesh.ChamberAlpha   = 0.3;
Mesh.ChamberRGB     = [0,0,1];


%================ Draw figure
LoadingFig          = EN_About(1);                                            	% Open loading screen
Fig.scnsize         = get(0,'ScreenSize');                                     	% Get screen resolution
Fig.Rect            = [0 0 Fig.scnsize(3), Fig.scnsize(4)];                    	% Set figure winow to fullscreen
Fig.FontSize        = 16;                                                     	% Set defualt font size
Fig.Background      = repmat(0.75,[1,3]);                                      	% Set figure background color  
Fig.InputBackground = repmat(0.85,[1,3]);                                       % Set input box background color
Fig.AxesBkgColor    = repmat(0.75,[1,3]);                                     	% Set axes background color
Fig.Handle          = figure('Name',sprintf('ElectroNav%c - Implant Vizualization',char(169)),... 	% Open a figure window with specified title
                    'Color',Fig.Background,...                                  % Set the figure window background color
                    'Renderer','OpenGL',...                                     % Use OpenGL renderer
                    'Position', Fig.Rect,...                                    % position figure window to fit fullscreen
                    'visible','off',...                                         % Figure remains invisible until complete
                    'NumberTitle','off',...                                     % Remove figure number from title
                    'IntegerHandle','off');                                     % Don't use integer handles
Fig.AxesH           = axes('position', [0.3, 0.1, 0.65, 0.8]);                  % Position main axes
Fig.CurrentMesh     = 1;                                                        

for i = 1:numel(Mesh.Filenames)
    if ~isempty(Mesh.Filenames{i})
        [Mesh.FV(i).vertices, Mesh.FV(i).faces]= stlread(Mesh.Filenames{i});
        if ApplyOffset(i)==1
            Mesh.FV(i).vertices = Mesh.FV(i).vertices+repmat(SL_offset,[size(Mesh.FV(i).vertices,1),1]);
        elseif ApplyOffset(i) == 2
            Mesh.FV(i).vertices = ENT_ApplyTform(ChamberTform, Mesh.FV(i).vertices);
            Mesh.FV(i).vertices = Mesh.FV(i).vertices+repmat(ChamberOffset,[size(Mesh.FV(i).vertices,1),1]);
        elseif ApplyOffset(i) == 3
            Mesh.FV(i).vertices = Mesh.FV(i).vertices*TowerTform;
            Mesh.FV(i).vertices = ENT_ApplyTform(ChamberTform, Mesh.FV(i).vertices);
            Mesh.FV(i).vertices = Mesh.FV(i).vertices+repmat(ChamberOffset,[size(Mesh.FV(i).vertices,1),1]);
        end
        Mesh.Handles{i} = patch('vertices',Mesh.FV(i).vertices, 'faces',Mesh.FV(i).faces,'edgecolor','none','facecolor',Mesh.FaceColor(i,:),'facealpha',Mesh.FaceAlpha(i));
    end
end

%=========== DRAW CYLINDER OF ACCESSIBLE TISSUE 
ChamberDiameter = 19;
[X,Y,Z]         = cylinder(ChamberDiameter/2, 100);
Z               = (Z*60)-30;
fvc             = surf2patch(X,Y,Z);
fvc.vertices    = ENT_ApplyTform(ChamberTform, fvc.vertices);
fvc.vertices    = fvc.vertices+repmat(ChamberOffset,[size(fvc.vertices,1),1]);
Mesh.ChamberH   = patch('Faces',fvc.faces,'Vertices',fvc.vertices,'FaceColor',Mesh.ChamberRGB,'facealpha',Mesh.ChamberAlpha, 'EdgeColor','none');
        
%=========== SET APPEARANCE
view(120,30);
axis equal;
grid on;
lighting phong
Fig.LighH = camlight('infinite');
set(gca,'fontsize',12,'tickdir','out','xtick',Fig.Xlims(1):10:Fig.Xlims(2),'ytick',Fig.Ylims(1):10:Fig.Ylims(2),'ztick',-30:10:90,'ylim',[-40, 60],'zlim',[-30,90]);
xlabel('Left - Right (mm)','fontsize',Fig.FontSize)
ylabel('Posterior - Anterior (mm)','fontsize',Fig.FontSize)
zlabel('Inferior - Superior (mm)','fontsize',Fig.FontSize)

%=========== PLOT REFERENCE PLANES
Fig.PlaneH(1) = patch([Fig.Xlims, Fig.Xlims([2,1])], [Fig.Ylims([1,1]), Fig.Ylims([2,2])], repmat(0,[1,4]), 'b','edgecolor','b');
Fig.PlaneH(2) = patch([Fig.Xlims, Fig.Xlims([2,1])], [Fig.Ylims([1,1]), Fig.Ylims([2,2])], repmat(SL_offset(3),[1,4]), 'g','edgecolor','g');
set(Fig.PlaneH,'facealpha',Mesh.PlaneAlpha, 'visible', Mesh.PlaneOn);



%% ========================== ADD GUI ELEMENTS ============================
Fig.GUIPannel   = uipanel('BackgroundColor',Fig.Background,'Units','normalized','Position',[0.05,0.1,0.2,0.8]);
Fig.PanelNames  = {'Implant geometry', 'Appearance', 'Options'};                            % Set GUI pannel titles
BoxPos(1,:)     = [400, 600, 260, 240];                                                     % Pannel 1 = cell info
BoxPos(2,:)     = [400, 380, 260, 210];                                                     % Pannel 2 = atlas structures
BoxPos(3,:)     = [400, 20, 260, 350];                                                      % Pannel 3 = MRI
for i = 1:numel(Fig.PanelNames)
    Fig.Handles.UIpannel(i) = uipanel('Title',Fig.PanelNames{i},'FontSize',18,'BackgroundColor',Fig.Background,'Units','pixels','Position',BoxPos(i,:),'Parent',Fig.GUIPannel);
end

%=========== Implant geometry
Fig.Struct.ButtonDim    = [100 25];
Fig.Struct.LabelStrings = {'Chamber diameter','Color','Opacity','Smoothing','',''};
Fig.Struct.InputType    = {'popupmenu','PushButton','slider','slider', 'checkbox','checkbox'};
% Fig.Struct.InputStrings = {Structures.Names, [], [], [], 'Visible','Wireframe'};
% Fig.Struct.InputValue   = {Structures.CurrentStructure, [], Structures.Opacity(Structures.CurrentStructure), Structures.Smoothing(Structures.CurrentStructure), Structures.On(Structures.CurrentStructure), Structures.Wire(Structures.CurrentStructure)};
% Fig.Struct.ButtonPos    = [repmat(10,[numel(Fig.Struct.LabelStrings),1]), [0:30:((numel(Fig.Struct.LabelStrings)-1)*30)]'+10];
% Fig.Struct.ButtonPos    = Fig.Struct.ButtonPos(end:-1:1,:);
% 
% %=========== Appearance
% Fig.Struct.LabelStrings = {'Selected surface','Load','Color','Opacity','Visible','Surface'};
% Fig.Struct.InputType    = {'popupmenu','PushButton','PushButton','edit','checkbox','popupmenu'};
% Fig.Struct.InputStrings = {{Mesh.FileTypes}, [], [], Mesh.FaceAlpha(Fig.CurrentMesh), 'Visible','Wireframe'};
% Fig.Struct.InputValues  = {Fig.CurrentMesh, 0, }


close(LoadingFig);
set(Fig.Handle, 'visible', 'on');
drawnow


end

%% ========================= GUI CALLBACKS ================================


