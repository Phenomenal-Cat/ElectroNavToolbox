
%============================= SlicePlot3D.m ==============================
%
%



Method = 3;
if Method == 1
    PulvMesh = '/Volumes/APM_B/PulvinarReview/Images/3D reconstructions/LeftPulvinar_2Xscale.stl';
    [v, f, n, c, stltitle] = stlread(PulvMesh);
    FV.faces = f';
    FV.vertices = v';
elseif Method == 2
    PulvMesh = '/Volumes/projects/murphya/EN_Data/Atlases/inia19/VTKs/pulvinar.vtk';
    [v,f] = read_vtk(PulvMesh);
    FV.faces = f';
    FV.vertices = v';
elseif Method == 3
    PulvVol = '/Volumes/projects/murphya/EN_Data/Atlases/inia19/Structures/pulvinar.nii';
    PulvVol = load_nii(PulvVol);
%     PulvVol.img = PulvVol.img(1:PulvVol.hdr.hist.originator(1),:,:);      	% Remove one hemipshere
%     Vol = smooth3(PulvVol.img, 'gaussian', 5);      % Smooth volume
    Vol =  PulvVol.img;
    isoval = 2;
    FV = isosurface(Vol,isoval);                                            % Threshold isosurface
%     FV2 = isocaps(1:80, 1:size(Vol,2), 1:size(Vol,3), Vol, isoval);             	% Get isocaps
end

% % FV = smoothpatch(FV);
% MeshColor       = [0.5 0.5 1];
% MeshAlpha       = 0.4;
% MeshEdgeColor   = 'none';
% CapColor        = 'b';
% CapAlpha        = 0.8;
% CapEdgeColor   	= [0 0 1];
% 
% % ph1 = contourslice([],[],[], Vol, [80:2:94], [], [], [10,10]);
% % set(ph1,'EdgeColor',[0 0 1],'linewidth',2, 'facecolor',[0 0 1]);
% 
% hcap = patch(FV2,'FaceColor',CapColor,'EdgeColor',CapEdgeColor,'Facealpha',CapAlpha);
% hold on;
% % ph2 = patch('vertices',FV.vertices, 'faces',FV.faces,'facecolor',MeshColor,'edgecolor',MeshEdgeColor,'facealpha',MeshAlpha);
% % set(ph2,'FaceVertexCData',repmat(MeshColor, size(FV.vertices)));
% % shading('interp');
% axis equal tight
% grid on
% % set(gca,'xticklabels',[],'yticklabels',[],'zticklabels',[]);
% set(gca,'tickdir','out');
% view(50,25);
% lh(1) = camlight('infinite');
% material([0.6,0.6,0.3,4,1])
% 
% set(gca,'xlim', [78, 94], 'ylim', [54, 82], 'zlim', [48 76]);

% lightangle(45,30); 
% lighting phong
% isonormals(Vol, ph2);


%================== ALTERNATIVE METHOD
figure;
SlicePos                = [80:2:94];
ThreshVol               = Vol;
MaxAlpha                = 0.6;
ThreshVol(ThreshVol<5)  = 0;
ThreshVol(ThreshVol>=5) = 1;
for s = 1:numel(SlicePos)
    Im = rot90(squeeze(ThreshVol(:,SlicePos(s),:)),3);
    ImH(s) = surf(repmat(SlicePos(s),[2,2]), [0 1; 0 1], [0 0; 1 1], 'FaceColor','texturemap','cdata', Im, 'FaceAlpha','texturemap','AlphaData',Im*MaxAlpha, 'AlphaDataMapping','scaled');
    hold on;
end

