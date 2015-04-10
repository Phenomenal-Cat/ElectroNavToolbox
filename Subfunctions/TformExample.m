
% T-form Example
% http://blogs.mathworks.com/steve/2006/08/17/spatial-transformations-three-dimensional-rotation/


addpath(genpath(cd));
GridFile = 'Grid/Grid1.stl';
[v, f, n, c, stltitle] = stlread(GridFile);

set(0, 'DefaultFigureRenderer','zbuffer');

BrainACPCNii = '/PROJECTS-1/murphya/Toolboxes/ElectroNavToolbox/Subjects/Layla/Layla_GridScan_ACPC.nii';
GridNii = '/PROJECTS-1/murphya/Toolboxes/ElectroNavToolbox/Grid/Grid1.nii';


%============== Load transformation matrix
% XformFile = '/PROJECTS/murphya/Toolboxes/ElectroNavToolbox/Subjects/Layla/Layla_GridScan_ACPC.xform';
XformFile = '/PROJECTS-1/murphya/Toolboxes/ElectroNavToolbox/Subjects/Layla/Layla_ACPC_to_grid_2.xform';
fileID = fopen(XformFile);                            	% Open Xform file
Xform = cell2mat(textscan(fileID,'%f %f %f %f\n'));    	% Read Xform to matrix
InverseXform = inv(Xform);                              % Calculate inverse transform matrix
fclose(fileID);        


%============= Load nifti
nii = load_nii(BrainACPCNii);
Grid = load_nii(GridNii);

GridOrigin = Grid.hdr.hist.originator(1:3);
GridVolumeSize = size(Grid.img);
GridPixSize = Grid.hdr.dime.pixdim(2:4);

VolumeSform = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];
VolumeInvSform = inv(VolumeSform);
VolumeSize = size(nii.img);
Origin = nii.hdr.hist.originator;

GridVolume = zeros(size(nii.img));
Indx = [Origin(1)-GridOrigin(1), Origin(1)+(GridVolumeSize(1)-GridOrigin(1))-1;...
               Origin(2)-GridOrigin(2), Origin(2)+(GridVolumeSize(2)-GridOrigin(2))-1;...
               Origin(3)-GridOrigin(3), Origin(3)+(GridVolumeSize(3)-GridOrigin(3))-1];
GridVolume(Indx(1,1):Indx(1,2), Indx(2,1):Indx(2,2), Indx(3,1):Indx(3,2)) = Grid.img;


ScaleMat = zeros(4,4);
ScaleMat(1,1) = 0.25; ScaleMat(2,2) = 0.25; ScaleMat(3,3) = 0.25; ScaleMat(4,4) = 1; 
ScaleMat(1:3,4) = 4;

VolumeSform2 = VolumeSform;
VolumeSform2(1:3,4) = 0.25;
VolumeSform2(1,1) = 1; VolumeSform2(2,2) = 1; VolumeSform2(3,3) = 1;
% VolumeSform2(4,1:3) = 1;

% Xform(1:3,4) = Xform(1:3,4)./GridPixSize';  % Convert translation from mm to voxels!
T = Xform;% / VolumeSform2;        % Faster, more accurate equivalent of Xform * inv(VolumeSform)         
% T = VolumeSform \ Xform;      % Faster, more accurate equivalent of inv(VolumeSform) * Xform
% T(1:3,4) = T(1:3,4)*4;          
% T(3,4) = -T(3,4);


%=========== View grid on brain slice
MLlim = [-Origin(1), VolumeSize(1)-Origin(1)]*0.25;
APlim = [-Origin(2), VolumeSize(2)-Origin(2)]*0.25;
ISlim = [-Origin(3), VolumeSize(3)-Origin(3)]*0.25;


imagesc(APlim(1):APlim(2), MLlim(1):MLlim(2), nii.img(:,:,160))
hold on;
gh = imagesc(APlim(1):APlim(2), MLlim(1):MLlim(2), GridVolume(:,:,160));
alpha(gh,GridVolume(:,:,160))
daspect([1 1 1]);
plot(xlim, [0 0], '-r', [0 0], ylim, '-r');

% GridVolume = affine3d(GridVolume, inv(VolumeSform), 1:size(GridVolume,1), 1:size(GridVolume,2), 1:size(GridVolume,3), 'linear');

% GridVolume2 = affine3d(GridVolume, InverseXform, 1:size(GridVolume,1), 1:size(GridVolume,2), 1:size(GridVolume,3), 'linear');
% GridVolume3 = affine3d(GridVolume, Xform, 1:size(GridVolume,1), 1:size(GridVolume,2), 1:size(GridVolume,3), 'linear');



% thetaZ = pi/32;              % First rotate phi radians about the z axis
% thetaX = pi/4;              % Second rotate theta radians about the former x axis
% thetaY = 0;                 % Third rotate psi radians about the former  
% scale = 1;
% translation = [0 -40 40];   


% % Make 3D shape
% [x,y,z] = ndgrid(-1:.025:1);
% blob = z <= 0 & z >= -0.75 & x.^2 + y.^2 <= sqrt(0.25);
% blob = blob | (z > 0 & (abs(x) + abs(y) <= (0.5 - z)));

% Make affine tform struct
% blob_center = (size(blob) + 1) / 2;


%========= T-form option #1
% if exist('angle2dcm.m','file')==2
%     T = zeros(4,4);
%     T(1:3,1:3) = angle2dcm(thetaZ, thetaY, thetaX, 'ZYX');
%     T(:,4) = [translation, 1];
% else

%========= T-form option #2
%     T1 = makehgtform('zrotate',thetaZ);
%     T2 = makehgtform('yrotate',thetaY);
%     T3 = makehgtform('xrotate',thetaX);
%     T4 = makehgtform('translate',translation);
%     T = T1 * T2 * T3;
%     T(:,4) = T4(:,4);


%% =================== Display MRI slices
figure;
Xoffset = 46;       % Offset in voxels from ACPC origin to grid origin

%========== SAGITTAL SLICE
xPos = repmat(MLlim(1),[2,2]);
yPos = [APlim(1), APlim(1); APlim(2), APlim(2)];
zPos = [ISlim; ISlim];
SagSlice = surf(xPos,yPos,zPos,'CData',double(squeeze(nii.img(Origin(1)-Xoffset,:,:))),'FaceColor','texturemap','EdgeColor','none');
hold on;

%========= CORONAL SLICE
xPos = [MLlim(2), MLlim(2); MLlim(1), MLlim(1)];
yPos = repmat(APlim(1),[2,2]);
zPos = [ISlim; ISlim];
CorSlice = surf(xPos,yPos,zPos,'CData',double(squeeze(nii.img(:,Origin(2),:))),'FaceColor','texturemap','EdgeColor','none');

%========= AXIAL SLICE
xPos = [MLlim(1), MLlim(1); MLlim(2), MLlim(2)];
yPos = [APlim; APlim];
zPos = repmat(ISlim(1),[2,2]);
AxSlice = surf(xPos,yPos,zPos,'CData',double(squeeze(nii.img(:,:,Origin(3)+40))),'FaceColor','texturemap','EdgeColor','none');

xlabel('x/ M-L');
ylabel('y/ A-P');
zlabel('z/ I-S');


%======= Plot origin
set(gca, 'xlim', MLlim, 'ylim', APlim, 'zlim', ISlim);
plot3(xlim, [0 0], [0 0], '-r', [0 0], [0 0], zlim, '-r', [0 0], ylim, [0 0], '-r');
view(120, 30);



%% =================== Display original mesh

if size(v,2)==3
   v(:,4) = 1;   
end
VolumeSform(1:3,4) = 1;     % Remove translation components from Sform matrix?
% v = VolumeSform\v';


v = v';
p = patch('faces',f,'vertices',v(1:3,:)','facecolor','r','EdgeColor', 'none');
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1]);
camlight
lighting gouraud
grid on;



%============ Manually transformed mesh
Translation = [-11.5;-30.75;29]; 
theta = [36; 0; -1];                                                                             % Rotation about x axis (degrees)
% v2 = v(1:3,:)'* [1,0,0;0,cosd(theta(1)),sind(theta(1));0,-sind(theta(1)),cosd(theta(1))];       % 1) apply rotation
% v2 = v2'+repmat(Translation,[1,size(v2,1)]);                                                    % 2) apply translation
% p2 = patch('faces',f,'vertices',v2(1:3,:)','facecolor','y','EdgeColor', 'none');
% daspect([1 1 1]);
% view(3)
% camlight
% lighting gouraud




%=========== Matrix transformed mesh 1
RotationSequence = 'ZYX';
theta = -theta([3,2,1]);
theta = theta/360*2*pi;
dcm = angle2dcm2(theta(1),theta(2),theta(3),RotationSequence);
dcm(:,4) = Translation;
dcm(4,:) = ones(1,4);

v3 = dcm*v;       % NB: transform matrix T must come FIRST in the multiplication with the 4 x n coordinate array!!!
p2 = patch('faces',f,'vertices',v3(1:3,:)','facecolor','g','EdgeColor', 'none');
daspect([1 1 1]);
view(3)
camlight
lighting gouraud
view(150,30);




% %========= Check with inverse transform...
% InvT = inv(T);
% v3 = InvT*v;
% p3 = patch('faces',f,'vertices',v3(1:3,:)','facecolor','y','EdgeColor', 'none');


colormap gray
set(gca,'fontsize',16)
set(get(gca,'xlabel'),'fontsize',20)
set(get(gca,'ylabel'),'fontsize',20)
set(get(gca,'zlabel'),'fontsize',20)


