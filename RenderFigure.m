function [] = RenderFigure

% Pulvinar targeting figure

Target = [0 0];
TargetDepth = -50;
Electrode.CurrentDepth = -35;

Brain.Specular = 0.5;
Brain.Ambient = 0.2;
Brain.Diffuse = 0.6;
Brain.Alpha = 1;%0.4;
Brain.RGB = [0.5 0.5 0.5];
Brain.DefaultView = [-120 20];
% Brain.ChamberAngle = 38;
Brain.ChamberAngle = 0;

Surface.Atlas = 'NeuroMaps';

switch Surface.Atlas
    case 'Paxinos'
        Surface.VTKfile = 'Niftii/Frey/Paxinos_surface.vtk';
     	Surface.Xlim = [-32 32];
        Surface.Ylim = [-60 30];
        Surface.Zlim = [-20 40];
        Brain.ChamberOrigin = [-11,-16,38];
        
    case 'NeuroMaps'
        Surface.VTKfile = 'Atlases/NeuroMaps/NeuroMapsVTKs/Cortical_surface.vtk';
        Surface.PulvVTK = 'Atlases/NeuroMaps/NeuroMapsVTKs/pulvinar.vtk';
       	Surface.Xlim = [-32 32];
        Surface.Ylim = [-60 30];
        Surface.Zlim = [-25 40];
%         Brain.ChamberOrigin = [-11,-33,24];
        Brain.ChamberOrigin = [-11, -4, 26]; 
        
        l_medial_pulvinar_nucleus  = 92;
        l_inferior_pulvinar_nucleus  = 95;
        l_lateral_pulvinar_nucleus  = 96;
        
    case 'Saleem-Logothetis'
        Surface.VTKfile = 'Niftii/McLaren/McLaren_surface.vtk';
        Surface.Xlim = [-32 32];
        Surface.Ylim = [-30 60];
        Surface.Zlim = [-10 50];
        Brain.ChamberOrigin = [-11,-16,38];
end 

%=========== Load grid
Grid = ENT_GetGridParams('19mm_cylindrical');
for dim = 1:3
    GridFV.vertices(:,dim) = Grid.vertices(:,dim)+Brain.ChamberOrigin(dim);
end
GridFV.faces = Grid.faces;



Grid.HoleDiameter = 0.5;                                                    % Grid hole diameter (mm)
Grid.InterHoleSpacing = 1;                                                  % Distance between centres of adjacent holes (mm)
Grid.HolesPerColumn = [5 9 11 13 15 15 17 17 17 17 17 15 15 13 11 9 5];     % Number of holes per column
Grid.HolesPerDim = numel(Grid.HolesPerColumn);                              
Grid.TotalHoles = sum(Grid.HolesPerColumn);                                 % Total number of grid holes
Grid.OuterRadius = Grid.InterHoleSpacing*(Grid.HolesPerDim+1)/2;            % Outer radius of grid
Grid.Width = 10;  
Grid.Height = 10;
Guide.Top = 10;
Grid.RGB = [1 1 0];

Electrode.Type = 'NN32';
switch Electrode.Type
    
    case 'PX24'                    %============ Plexon V-probe 24 channel linear multi-electrode array (microfil)
     	Electrode.Length = 90;                  % Full electrode shaft length (mm)
        Electrode.Diameter = 0.2;               % shadt diameter (mm)
        Electrode.TipLength = 1.5;              % distance from tip to first contact (mm)
        Electrode.ContactSpacing = 0.2;         % distance between adjacent contacts (mm)
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization
        Electrode.ContactNumber = 24;           % total number of contacts
        Electrode.ContactLength = Electrode.ContactSpacing*Electrode.ContactNumber;
        Electrode.Colour = [1 1 0];             
        Electrode.ContactColour = [1 0 0];  
        
    case 'AO24'                     %============ Alpha Omega 24 channel linear multi-electrode array (microfil)
        Electrode.Length = 70;                  % Full electrode shaft length (mm)
        Electrode.Diameter = 0.3;               % shadt diameter (mm)
        Electrode.TipLength = 1.5;              % distance from tip to first contact (mm)
        Electrode.ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization
        Electrode.ContactNumber = 24;           % total number of contacts
        Electrode.ContactLength = Electrode.ContactSpacing*Electrode.ContactNumber;
        Electrode.Colour = [0 1 0];             
        Electrode.ContactColour = [1 0 0];      
        
    case 'NN32'                     %============ NeuroNexus 32 channel linear multi-electrode array
        Electrode.Length = 60;
        Electrode.Diameter = 0.3;
        Electrode.TipLength = 1;
        Electrode.ContactSpacing = 0.3;
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization
        Electrode.ContactNumber = 32;
        Electrode.ContactLength = Electrode.ContactSpacing*Electrode.ContactNumber;
        Electrode.Colour = [0 0 1];
        Electrode.ContactColour = [1 0 0];
end

[v,f] = read_vtk(Surface.VTKfile);
FV.vertices = v';
FV.faces = f';
FV.facevertexcdata = Brain.RGB;
FV.facecolor = 'flat';
FV.facealpha = Brain.Alpha;
FV.edgecolor = 'none';    
Brain.Object = patch(FV,'EdgeColor','none');
hold on;
% camlight above;
% camlight headlights;
lighting phong;
colormap bone;
axis(gca,'vis3d');                                      % Maintain axes ratio (do not scale)     
% Brain.Labels(1) = xlabel('Medial-Lateral','fontsize',16);                                        
% Brain.Labels(2) = ylabel('Posterior-Anterior','fontsize',16);
% Brain.Labels(3) = zlabel('Inferior-Superior','fontsize',16);
grid on;
axis equal;
view(Brain.DefaultView);
set(Brain.Object,'SpecularStrength',Brain.Specular,'AmbientStrength',Brain.Ambient,'DiffuseStrength',Brain.Diffuse);

if isfield(Surface,'PulvVTK')
    [v,f] = read_vtk(Surface.PulvVTK);
    FV.vertices = v';
    FV.faces = f';
    FV.facevertexcdata = [1 0 0];
    FV.facecolor = 'flat';
    FV.facealpha = 1;
    FV.edgecolor = 'none';    
    Pulv.Object = patch(FV,'EdgeColor','none');
end
set(Pulv.Object,'Facecolor',[1 0 0]);


%=============== Draw chamber & grid holes
Brain.Chamber = patch(GridFV, 'FaceColor', Grid.RGB, 'EdgeColor', 'none');
rotate(Brain.Chamber,[1,0,0],Brain.ChamberAngle, Brain.ChamberOrigin);
rotate(Brain.Chamber,[0,1,0],-10, Brain.ChamberOrigin);
set(Brain.Chamber,'FaceLighting','phong');%'FaceColor','interp',



[X,Y,Z] = cylinder(Electrode.Diameter/2,100);
X = X+Brain.ChamberOrigin(1)+Target(1);
Y = Y+Brain.ChamberOrigin(2)+Target(2);
Z1 = (Z*(Electrode.Length-Electrode.TipLength))+Electrode.CurrentDepth+Electrode.TipLength+Brain.ChamberOrigin(3);
E(1) = mesh(X,Y,Z1,Electrode.Colour,'FaceLighting','phong','EdgeColor','none');
set(E,'Facecolor',[0 0 1]);
rotate(E(1),[1,0,0],Brain.ChamberAngle, Brain.ChamberOrigin);


set(Brain.Object,'FaceVertexAlphaData',Brain.Alpha);
set(Brain.Object,'FaceAlpha',Brain.Alpha);
set(Brain.Object,'Facecolor',Brain.RGB);
whitebg([1 1 1]);
set(gca,'xlim',Surface.Xlim);
set(gca,'ylim',Surface.Ylim);
set(gca,'zlim',Surface.Zlim);
axis off;

% set(gcf,'Color',[0 0 0]);

% Settings
set(Brain.Object,'Visible','on');
set(Brain.Object,'HitTest','on');
set(gcf,'Renderer','opengl')
drawnow


% FileName = 'Pulvinar_render3';
% res_dpi = 300;                              
% set(gcf,'InvertHardcopy','off');
% print(gcf, sprintf('-r%d', res_dpi), '-dtiff', [FileName,'.tif']);



%================= Capture movie?
CamDur = 4;                                             % Animation duration (seconds)
CamStart = [45,20,1];                                	% Camera start position (azimuth, elevation, zoom) degrees
CamEnd = [405,20,1.02];                                	% Camera end position (azimuth, elevation, zoom) degrees
CamInc = (CamEnd-CamStart)/CamDur;                      % Camera movement increments (degrees per second)
            
% % for a = 1:360
% %     view(a,20);
% %     drawnow;
% % %     M(a) = getframe;
% % 
% % end
ViewZ = [CamStart([1,2]);CamEnd([1,2])];
OptionZ.FrameRate = 30;
OptionZ.Duration = 4;
FileName = 'ElectrodeMovie';
% CaptureFigVid(ViewZ, FileName,OptionZ);



%=============== Setup view for render
view(-90,0);
lh(1) = camlight('headlight');
lh(2) = camlight(0, 70);
grid on
axis on
set(gca,'fontsize',14);
set(gca,'ylim',[-60 40],'ytick',[-60:10:60],'zlim',[-40 40],'ztick',-40:10:40);

% Draw origin
sliceh = plot3(repmat(-30,[1,2]),[30,-50],[0, 0],'-g','linewidth',1);
sliceh = plot3(repmat(-30,[1,2]),[0,0],[-30, 30],'-g','linewidth',1);

% Draw slice
PulvSlice = -12;
STSslice = -4;
sliceh = plot3(repmat(-30,[1,2]),[STSslice,STSslice],[-30, 30],'-r','linewidth',3);
axis off

Filename = '/Users/aidanmurphy/Pulvinar_project/ChunshansRenders/STSLateralView3D.png';
export_fig(Filename, '-png','-transparent','-m3');



end


%% ========================= SUBFUNCTIONS =================================
function h = PlotCircle(x,y,z,r,c)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
zunit = repmat(z,size(th));
h = plot3(xunit, yunit,zunit,c);
end

function h = FillCircle(target,z,r,N,c)
THETA=linspace(0,2*pi,N);
RHO=ones(1,N)*r;
[X,Y] = pol2cart(THETA,RHO);
X=X+target(1);
Y=Y+target(2);
Z = repmat(z,[1,N]);
h=fill3(X,Y,Z,c,'EdgeColor','none');
end

function h = PlotSphere(x,y,z,r,c)
[X,Y,Z] = sphere(100);
X = (X*r)+x;
Y = (Y*r)+y;
Z = (Z*r)+z;
h=mesh(X,Y,Z,'FaceColor',c,'EdgeColor','none');

% ellipsoid(x,y,z,1,1,Electrode.ContactLength)

end

%============================== DRAW GRID =================================
function GridObject = DrawGrid(Grid,Position)
    GridObject(1) = FillCircle(Position([1,2]),Position(3),Grid.OuterRadius,100,'y');
    hold on;
    GridObject(end+1) = FillCircle(Position([1,2]),Position(3)+Grid.Height,Grid.OuterRadius,100,'y');
  	[X,Y,Z] = cylinder(Grid.OuterRadius,1000);
    X = X+Position(1);
 	Y = Y+Position(2);
    Z = (Z*Grid.Height)+Position(3);
    GridObject(end+1) = mesh(X,Y,Z,'FaceColor',Grid.RGB,'EdgeColor','none');
    for z = [Grid.Height+Position(3), Position(3)]
        for i = 1:Grid.HolesPerDim
            for h = 1:Grid.HolesPerColumn(i)
                x = ((-((Grid.HolesPerColumn(i)-1)/2)+(h-1))*Grid.InterHoleSpacing)+Position(1);
                y = ((((Grid.HolesPerDim+1)/2)-i)*Grid.InterHoleSpacing)+Position(2);
                r = Grid.HoleDiameter/2;
                GridObject(end+1) = PlotCircle(x,y,z,r,'k');
            end
        end
        GridObject(end+1) = plot3([Position(1) Position(1)],[-Grid.OuterRadius,Grid.OuterRadius]+Position(2),[z z],'-k');
        GridObject(end+1) = plot3([-Grid.OuterRadius,Grid.OuterRadius]+Position(1),[Position(2) Position(2)],[z z],'-k');
    end
    GridParent = hgtransform('Parent',gca);
    set(GridObject,'Parent',GridParent);
    
    
end

