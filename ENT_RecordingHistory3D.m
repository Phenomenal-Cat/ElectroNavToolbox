function fh = ENT_RecordingHistory3D(SubjectID, Structures)

if nargin == 0
    SubjectID = 'Dexter';
end
Radius      = 0.1;

[ContactCoords, SessionParams] = EN_GetContactCoordinates([], SubjectID);

%========== Plot all electrode locations
fh = figure;
Eh = [];
PenetrationColors = jet(size(ContactCoords,1));
for d = 1:size(ContactCoords,1)
    for ch = 1:size(ContactCoords,3)
        Eh(end+1) = PlotSphere(squeeze(ContactCoords(d,:,ch)), Radius, PenetrationColors(d,:));
        hold on;
    end
end

%========== Plot reference structures
if ~exist('Structures', 'var')
    Defaults        = ENT_LoadDefaults(SubjectID);
    AllStructures   = wildcardsearch(Defaults.VTKdir, '*.vtk');
    [Selection,ok]  = listdlg('ListString',AllStructures,'SelectionMode','multi','PromptString','Select structures:','ListSize',[400, 200]);
    Structures      = AllStructures(Selection);
end
MeshH               = ENT_Plot3DMeshes(AllStructures(Selection));


%========== Format plot
grid on;
lighting phong
camlight('infinite');
daspect([1 1 1]);
xlabel('M-L (mm)')
ylabel('P-A (mm)')
zlabel('I-S (mm)')

cbh = colorbar;
colormap jet;
set(cbh, 'ticks', linspace(-1, 1, numel(SessionParams)), 'ticklabels', {SessionParams.DateString});

set(gca, 'xlim',[0 18], 'ylim', [-20 -4], 'zlim', [-8 12]);
title(sprintf('%s recording history %s - %s', SubjectID, SessionParams(1).DateString, SessionParams(end).DateString), 'fontsize', 18);

end


function h = PlotSphere(XYZ, r, c)

N       = 20;
[X,Y,Z] = ellipsoid(XYZ(1),XYZ(2), XYZ(3), r, r, r, N);
h       = surface(X,Y,Z,repmat(0, size(Z)),'facecolor',c,...
                'EdgeColor','none','FaceAlpha',0.8,'ButtonDownFcn',@ContactSelect); % Apply value to render    
end