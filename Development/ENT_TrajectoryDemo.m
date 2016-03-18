function EN_TrajectoryDemo(SubjectID)

%========================= EN_TrajectoryDemo.m ============================
% This function renders and saves an animation to illustrate the trajectory
% of the recording grid relative to target structures, and the electrode
% positions across sessions.
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy
% ? Copyleft 2014, GNU General Public License
%==========================================================================

global ChamberTform

if nargin == 0
    SubjectID   = 'Layla';
end

%================== LOAD DATA
Defaults    = ENT_LoadDefaults(SubjectID);          % Get default files for subject
load(Defaults.Xform);                               % Load transform matrix
SessionParams = ENT_LoadSessionParams(Defaults.HistoryFile);
GridOffset 	= T*[0 0 0 1]';                         % get grid origin offset from AC

%================== LOAD MRI
Nii         = load_nii(Defaults.MRI);
Yrange      = ([0, Nii.hdr.dime.dim(3)]-Nii.hdr.hist.originator(2))*Nii.hdr.dime.pixdim(3);
Zrange      = ([0, Nii.hdr.dime.dim(4)]-Nii.hdr.hist.originator(3))*Nii.hdr.dime.pixdim(4);
SliceIndx   = Nii.hdr.hist.originator(1)+round(GridOffset(1)/Nii.hdr.dime.pixdim(2));
Slice       = squeeze(Nii.img(SliceIndx,:,:));
Slice       = double(Slice);
Slice(Slice>15000) = 15000;
Slice       = Slice/15000;

switch SubjectID
    case 'Layla'
        SliceXlim = [-16, 0];
        SliceZlim = [-10, 15];
        MRIPlanePos = -30;
    case 'Dexter'
        SliceXlim = [0, 16];
        SliceZlim = [-10, 15];
        MRIPlanePos = 0;
end


MeshFiles{1} = '/Volumes/PROJECTS/murphya/EN_data/Subjects/Dexter/VTKs/pulvinar.vtk';
MeshFiles{2} = '/Volumes/PROJECTS/murphya/EN_data/Subjects/Dexter/VTKs/Cortical_surface.vtk';
% MeshFiles{3} = '/Volumes/APM_B/ANATOMY/SurgeryPlanning/SurfaceRecon/Dexter_scalp_surf2.stl';
MeshFiles{4} = '/Volumes/PROJECTS/murphya/ANATOMY/SurgeryPlanning/SurfaceRecon/Grid1.stl';
MeshFiles{3} = '/Volumes/PROJECTS/murphya/ANATOMY/SurgeryPlanning/SurfaceRecon/Grid1.stl';
% MeshFiles{5} = '/Volumes/APM_B/ANATOMY/SurgeryPlanning/SurfaceRecon/FlexMT_scaled2.stl';


ApplyOffset     = [1,1,1,2,3];      
TowerTform      = makehgtform('zrotate',pi/2);      
TowerTform      = TowerTform(1:3,1:3);
ChamberTform   	= T;
FaceAlpha       = [1, 0.5, 0.1, 0.5, 0.25];
FaceColor       = [1 0 0; 1 0.7 0.7; 0.5 0.5 0.5; 1 1 0; 1 1 0];

%================== PLOT MAIN FIGURE
Fig.Rect    = get(0,'ScreenSize');                            % Specify full screen square
Fig.Rect    = [-1919,1,1920,1007];
Fig.fh      = figure('position',Fig.Rect, 'color', [1 1 1]);
Fig.axh(1)  = axes('position',[0.3 0.1 0.7 0.85]);
for i = 1:numel(MeshFiles)
    switch MeshFiles{i}(end-2:end)
        case 'stl'
            [FV(i).vertices, FV(i).faces]= stlread(MeshFiles{i});
        case 'vtk'
            [FV(i).vertices, FV(i).faces]= read_vtk(MeshFiles{i});
            FV(i).vertices = FV(i).vertices';
            FV(i).faces = FV(i).faces';
    end
  	if ApplyOffset(i) == 2
        FV(i).vertices = ENT_ApplyTform(ChamberTform, FV(i).vertices);
    elseif ApplyOffset(i) == 3
        FV(i).vertices = FV(i).vertices*TowerTform;
        FV(i).vertices = ENT_ApplyTform(ChamberTform, FV(i).vertices);
    end
    Mesh.Handles{i} = patch('vertices',FV(i).vertices, 'faces',FV(i).faces,'edgecolor','none','facecolor',FaceColor(i,:),'facealpha',FaceAlpha(i), 'edgealpha', FaceAlpha(i));
    hold on;
end
set(Mesh.Handles{2}, 'BackFaceLighting','reverselit')
% set(Mesh.Handles{1}, 'edgecolor', FaceColor(i,:))

%=========== DRAW MRI SLICE
PlaneH = surface(repmat(MRIPlanePos,[2,2]),[Yrange',Yrange'],[Zrange; Zrange]);
set(PlaneH,'edgeColor','none','facecolor','texturemap','cdata', repmat(Slice,[1,1,3]));
colormap gray

%=========== DRAW CYLINDER OF ACCESSIBLE TISSUE 
[X,Y,Z] = cylinder(8,100);
Z = (Z*80)-40;
fvc = surf2patch(X,Y,Z);
fvc.vertices = ENT_ApplyTform(ChamberTform, fvc.vertices);
% fvc.vertices = fvc.vertices+repmat(ChamberOffset,[size(fvc.vertices,1),1]);
GT = patch('Faces',fvc.faces,'Vertices',fvc.vertices,'FaceColor',[0 1 1],'facealpha',0, 'EdgeColor','none');
        

%=========== DRAW ALL ELECTRODES
for S = 1:numel(SessionParams)
    for n = 1:numel(SessionParams(S).Target)
        elh{S,n} = PlotElectrode(SessionParams(S).Target{n}, SessionParams(S).Depth{n});
        set(elh{S,n},'visible','off');
        hold on
    end
end

%=========== DRAW CORONAL SLICES
SlicePos = -18:1:-8;
for p = 1:numel(SlicePos)
    CorSlice(p) = surface([SliceXlim',SliceXlim'],repmat(SlicePos(p),[2,2]),[SliceZlim; SliceZlim]);
end
set(CorSlice,'edgeColor',[0 0 1],'facecolor',[0 0 1],'facealpha',0,'edgealpha',0);



%=========== SET APPEARANCE
view(90,0);
axis equal tight;
axis vis3d
grid on;
daspect([1,1,1]);
lighting phong
LH = camlight('infinite');
if numel(MeshFiles) >= 6
    set(gca,'fontsize',12,'tickdir','out','xtick',-40:10:40,'ytick',-60:10:70,'ztick',-30:10:200,'zlim',[-30,180]);
elseif numel(MeshFiles) <= 5
    set(gca,'fontsize',12,'tickdir','out','xtick',-40:10:40,'ytick',-50:10:50,'ztick',-30:10:90,'ylim',[-50, 40],'zlim',[-30,60]);
end
xlabel('Left - Right (mm)','fontsize',16)
ylabel('Posterior - Anterior (mm)','fontsize',16)
zlabel('Inferior - Superior (mm)','fontsize',16)


%=========== PLOT 2D GRID VIEW
Fig.axh(2)  = axes('position',[0.15, 0.3,0.2,0.4]);
[Grid.vertices, Grid.faces]= stlread(MeshFiles{4});
GridH = patch('vertices',Grid.vertices, 'faces',Grid.faces,'edgecolor','none','facecolor',[1 1 0],'facealpha',1);
axis equal;
grid on;
daspect([1,1,1]);
lighting phong
camlight('infinite');
set(Fig.axh(2), 'visible','off','xlim',[-10, 10], 'ylim',[-10,10]);
set(GridH, 'visible', 'off');
xlabel('Medial-lateral (mm)','fontsize',16);
ylabel('Posterior-anterior (mm)','fontsize',16);
hold on


% %=========== PLOT REFERENCE PLANES
% Xlims = get(gca,'xlim');
% Ylims = get(gca,'ylim');
% PlaneH(1) = patch([Xlims, Xlims([2,1])], [Ylims([1,1]), Ylims([2,2])], repmat(0,[1,4]), 'b','edgecolor','b');
% PlaneH(2) = patch([Xlims, Xlims([2,1])], [Ylims([1,1]), Ylims([2,2])], repmat(SL_offset(3),[1,4]), 'g','edgecolor','g');
% set(PlaneH,'facealpha',0.1);


set(Mesh.Handles{3}, 'visible', 'off');

%=========== PREPARE MOVIE CAPTURE
Filename            = sprintf('%s_ElectrodeTrajectoryDemo.mov', SubjectID);
CaptureRect         = get(gca,'position');
Compression         = 'Photo PNG';
Transparency        = true;
movObj              = QTWriter(Filename, 'MovieFormat', Compression, 'Transparency', Transparency);
movObj.FrameRate    = 60;  

%=========== RUN ANIMATION LOOP AND CAPTURE FRAMES
TotalDuration   = 4.5;
FrameRate       = 60;
NoFrames        = FrameRate*TotalDuration;
StartPenetrations = 2*FrameRate;
PenDur          = round((NoFrames-StartPenetrations)/numel(SessionParams))
Alpha.Plane     = ones(1, NoFrames);
Alpha.Grid      = zeros(1, NoFrames);
Alpha.Cortex    = zeros(1, NoFrames);
Alpha.Pulv      = zeros(1, NoFrames);
Alpha.GT        = zeros(1, NoFrames);
Background      = zeros(1, NoFrames);

Alpha.Plane(FrameRate*1+1:FrameRate*2) = linspace(1,0,FrameRate*1);
Alpha.Plane(FrameRate*2+1:end) = 0;
Background = ones(1, NoFrames)-Alpha.Plane;
Alpha.Grid(1:(FrameRate*1)) 	= linspace(0, 0.5, FrameRate*1);
Alpha.Grid(FrameRate*1+1:end)	= 0.5;
Alpha.Cortex(((FrameRate*1)+1):FrameRate*2) = linspace(0, 0.5, FrameRate*1);
Alpha.Cortex(FrameRate*2+1:end) = 0.5;
Alpha.Cortex(FrameRate*2.5+1:FrameRate*3.5) = linspace(0.5, 0, FrameRate*1);
Alpha.Cortex(FrameRate*3.5+1:end) = 0;
Alpha.Pulv(((FrameRate*1)+1):FrameRate*2) = linspace(0, 0.8, FrameRate*1);
Alpha.Pulv(FrameRate*2+1:end) = 0.8;
Alpha.Pulv(FrameRate*2.5+1:end) = linspace(0.8, 0.2, FrameRate*2);


%=========== Schedule plot rotation and zoom
TargetAzEl      = [145, 15];
ZoomFramesIndx  = FrameRate*1.5+1:FrameRate*TotalDuration;
Views.Az        = repmat(90, [1, NoFrames]);
Views.El        = zeros([1, NoFrames]);
Views.Az(ZoomFramesIndx) = linspace(90, TargetAzEl(1), numel(ZoomFramesIndx));
Views.El(ZoomFramesIndx) = linspace(0, TargetAzEl(2), numel(ZoomFramesIndx));
Views.Az(FrameRate*4+1:end) = linspace(Views.Az(FrameRate*4), 90, FrameRate*0.5);
Views.El(FrameRate*4+1:end) = linspace(Views.El(FrameRate*4), 0, FrameRate*0.5);
AxLimTargets = [SliceXlim; -18 -8; -10 15];
AxLims{1} = repmat(get(Fig.axh(1),'xlim'), [NoFrames, 1]);
AxLims{2} = repmat(get(Fig.axh(1),'ylim'), [NoFrames, 1]);
AxLims{3} = repmat(get(Fig.axh(1),'zlim'), [NoFrames, 1]);
for a = 1:numel(AxLims)
    AxLims{a}(ZoomFramesIndx,1) = linspace(AxLims{a}(1,1), AxLimTargets(a,1), numel(ZoomFramesIndx));
    AxLims{a}(ZoomFramesIndx,2) = linspace(AxLims{a}(1,2), AxLimTargets(a,2), numel(ZoomFramesIndx));
end
CorSliceFaceAlpha = zeros(1, NoFrames);
CorSliceEdgeAlpha = zeros(1, NoFrames);
CorSliceFaceAlpha(4*FrameRate+1:end) = linspace(0,0.1,FrameRate*0.5);
CorSliceEdgeAlpha(4*FrameRate+1:end) = linspace(0,1,FrameRate*0.5);

S = 1;
for f = 1:NoFrames
    fprintf('Frame %d/%d (%.1f%%)....\n', f, NoFrames, (f/NoFrames)*100);
    set(PlaneH, 'facealpha', Alpha.Plane(f));
    set(Mesh.Handles{4}, 'facealpha', Alpha.Grid(f));
    set(Mesh.Handles{2}, 'facealpha', Alpha.Cortex(f));
    set(Mesh.Handles{1}, 'facealpha', Alpha.Pulv(f));
    set(GT, 'facealpha', Alpha.GT(f));
    set(CorSlice,'facealpha',CorSliceFaceAlpha(f),'edgealpha',CorSliceEdgeAlpha(f));
    axes(Fig.axh(1));                                                               
    view(Views.Az(f), Views.El(f));                                                 % Adjust camera angle
    set(gca, 'color', repmat(Background(f), [1,3]));                                % Adjust background color
    set(gca,'xlim',AxLims{1}(f,:),'ylim', AxLims{2}(f,:),'zlim',AxLims{3}(f,:));    % Adjust axes limits
    
    if f == FrameRate*2
        set(Fig.axh(2), 'visible','on');
        set(GridH, 'visible', 'on');
    end
    
    if f >= StartPenetrations+(PenDur*(S-1))
        if S < numel(SessionParams)
            for n = 1:numel(SessionParams(S).Target)
                axes(Fig.axh(2));
                ph(S,n) = plot3(SessionParams(S).Target{n}(1), SessionParams(S).Target{n}(2), 10, '.b','markersize', 40);
                set(elh{S,n},'visible','on');
            end
            title(sprintf('%S - Depth = %.2fmm', SessionParams(S).Date, SessionParams(S).Depth{1}),'fontsize', 18);
            S = S+1;
        end
    end
    drawnow
%     Frame = getframe(Fig.fh, Fig.Rect+[100 10 -200 -220]);   
    Frame = screencapture(Fig.fh); 
    Frame = Frame(80:end,100:end-100,:);
    writeMovie(movObj, Frame);  
    
end
close(movObj);


end



%============== PLOT 3D ELECTRODE 
function E = PlotElectrode(Target, Depth)
    global ChamberTform
    
    Electrode = ENT_GetElectrodeParams('PLX24');
    [X,Y,Z] = cylinder(Electrode.Diameter/2,100);           % Electrode shaft
    X = X+Target(1);
    Y = Y+Target(2);
    Z1 = (Z*(Electrode.Length-Electrode.TipLength))-Depth+Electrode.TipLength;
    fvc1 = surf2patch(X,Y,Z1);
    fvc1.vertices = ENT_ApplyTform(ChamberTform,fvc1.vertices);
    
    [X2,Y2,Z2] = cylinder([0 Electrode.Diameter/2]);        % Electrode tip
    X2 = X2+Target(1);
    Y2 = Y2+Target(2);
    Z2 = (Z2*Electrode.TipLength)-Depth;
    fvc2 = surf2patch(X2,Y2,Z2);
    fvc2.vertices = ENT_ApplyTform(ChamberTform,fvc2.vertices);
    
    [X3,Y3,Z3] = cylinder(Electrode.Diameter*0.55,100);     % Electrode contacts
    X3 = X3+Target(1);
    Y3 = Y3+Target(2);

    Electrode.E(1) = patch('Faces',fvc1.faces,'Vertices',fvc1.vertices,'FaceColor',Electrode.MRIColor,'EdgeColor','none');
   	Electrode.E(2) = patch('Faces',fvc2.faces,'Vertices',fvc2.vertices,'FaceColor',Electrode.MRIColor,'EdgeColor','none');
    NoCurrentContacts = 0;
    for c = 1:Electrode.ContactNumber                                       % Loop through max number
        Z3 = (Z*Electrode.ContactDiameter)-Depth+Electrode.TipLength+(c-1)*Electrode.ContactSpacing;
        [xa ya za] = ENT_ApplyTform(ChamberTform, X3(1,:),Y3(1,:),Z3(1,:));
        [xb yb zb] = ENT_ApplyTform(ChamberTform, X3(2,:),Y3(2,:),Z3(2,:));
        if c <= NoCurrentContacts && c <= Electrode.ContactNumber            % 1) Contact # already exists: move it
            set(Electrode.E(2+c), 'xdata', [xa; xb], 'ydata', [ya;yb], 'zdata', [za;zb]);
        elseif c > NoCurrentContacts && c <= Electrode.ContactNumber         % 2) Contact # doesn't exist: create it
            Electrode.E(2+c) = mesh([xa; xb], [ya;yb], [za;zb],'FaceColor',Electrode.ContactColor,'EdgeColor','none');
            hold on;
        elseif c > Electrode.ContactNumber && c <= NoCurrentContacts         % Contact numbers > requested number: delete it
            delete(Electrode.E(2+c));
        end
    end
    E = Electrode.E;
end
