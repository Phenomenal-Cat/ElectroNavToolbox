
%=========================== RenderMacaqueBrain3D.m =======================
% Render 3D surface of macaque template brain with slices and structures plotted
%
%
%
%==========================================================================


VTKdir = '/Volumes/APM_1/ElectroNavToolbox/Subjects/Layla/VTKs/';
VTKfiles = dir([VTKdir,'*.vtk']);

ENRootDir = '/Users/aidanmurphy/Pulvinar_project/ElectroNavToolbox';
Atlasfile = fullfile(ENRootDir,'/Atlases/NeuroMaps/inia19-NeuroMaps.nii');
MRIfile = fullfile(ENRootDir,'/Atlases/NeuroMaps/inia19-t1-brain.nii');
Skullfile = '/Volumes/APM_1/MacaqueAtlasTools/Niftii/Skulls/Skull_ITK_2.vtk';

Save = 0;                                   % Automatically save figure image?
smooth = 5;                                 % Smoothing kernel diameter (voxels)
thresh = 50;                                % -Inf = Otsu's method; Inf = median intensity
reduce = 0.15;                              % Reduce mesh complexity to X % of original      
    
FigureBackground = [0 0 0];
VTK.fc = {'w','g','b','m','c','r','y'};   	% Set surface mesh face color
VTK.ec = 'none';                          	% Set surface mesh edge color
VTK.bf = 'reverselit';                      % Set surface mesh backface lighting
VTK.alpha =  ones(1,numel(VTKfiles))*0.1;	% Set surface mesh face alpha


SliceView = -1;                             % <1 = 3D; 0 = axial; 1 = sagittal; 2 = coronal;
SliceColor = 'b';                           % Set slice edge color
SliceAlpha = 0.4;                           % Set slice alpha value
SlicePos = -13:1:-4;                        % Set slice positions (mm)


%============= General material lighting properties
Ambient = 0.3;                          
Diffuse = 0.5;
Specular = 0.4;
SpecExp = 6;
SpecCol = 1;
 


% %% ========================== READ MRI DATA ===============================
% if exist('spm_vol','file')==2               % If SPM functions are detected...
%     Hdr = spm_vol(MRIfile);               	% use SPM calls
%     Vol = spm_read_vols(Hdr);               
% else
%     nii = read_nii(MRIfile);                % Use niftii toolbox calls
%     Vol = nii.img;                          
%     
% end
% Vol(isnan(Vol)) = 0;                        % Convert any NaN vlaues to zero
% if (round(smooth) > 3)                    	% blur image prior to edge extraction
%     fprintf('Applying gaussian smooth with %d voxel diameter\n',round(smooth));
%     Vol = smooth3(Vol,'gaussian',round(smooth));
% end;
% if (isinf(thresh) && (thresh < 0))          % if thresh = -Inf, use Otsu's method
%      thresh = graythresh(Vol);
% %      thresh = otsu(Vol);                    % use Otsu's method to detect isosurface threshold
% elseif (isnan(thresh)) || (isinf(thresh))   % if +Inf, use midpoint
% 	thresh = max(Vol(:)) /2;                % use  max/min midpoint as isosurface threshold
% %     thresh = mean(Vol(:));                % use mean-sd to detect isosurface threshold - heavily influenced by proportion of dark air
% end;
% 
% FV = isosurface(Vol,thresh);
% if (reduce ~= 1.0)                          % next: simplify mesh
%     FV = reducepatch(FV,reduce);
% end;
% layer = 1;
% 
% % for layer = 1
%     v.surface(layer).faces = FV.faces;
%     v.surface(layer).vertices = FV.vertices;
%     v.vprefs.demoObjects = false;
%     clear('FV');
%     %next: isosurface swaps the X and Y dimensions! size(Vol)
%     i = 1;
%     j = 2;
%     v.surface(layer).vertices =  v.surface(layer).vertices(:,[1:i-1,j,i+1:j-1,i,j+1:end]);
%     %BELOW: SLOW for loop for converting from slice indices to mm
%     %for vx = 1:size( v.surface(layer).vertices,1) %slow - must be a way to do this with bsxfun
%     % wc = Hdr.mat * [ v.surface(layer).vertices(vx,:) 1]'; %convert voxel to world coordinates
%     % v.surface(layer).vertices(vx,:) = wc(1:3)';
%     %end
%     %BELOW: FAST vector for converting from slice indices to mm
%     vx = [ v.surface(layer).vertices ones(size( v.surface(layer).vertices,1),1)];
%     vx = mtimes(Hdr.mat,vx')';
%     v.surface(layer).vertices = vx(:,1:3);
% % end
% 
% 
% 
% %     Skull.v = verts';
% %     Skull.f = faces';
% % 
% %     CullIndx = find(verts < -700);
% %     Faces(CullIndx, :) = [];
% 
% % rotate(skullh,[1 0 0],-20);
% % rotate(skullh,[0 0 1], -90);



%% ================== Make new VTKs from atlas volume
AtlasNii = load_nii(Atlasfile);




%% ================== Load VTK surfaces
for v = 1:numel(VTKfiles)
    [verts,faces] = read_vtk(fullfile(VTKdir, VTKfiles(v).name));
    VTK.handle(v) = patch('vertices', verts','faces', faces',...
                            'edgecolor',VTK.ec,'BackFaceLighting',VTK.bf,...
                            'facealpha',VTK.alpha(v),'facecolor',VTK.fc{v},'facelighting','phong');
                        
 	if ~isempty(strfind(VTKfiles(v).name,'pulvinar'))
        set(VTK.handle(v), 'facealpha', 1);
    end
    hold on;
end
set(gca,'DataAspectRatio',[1 1 1]);
set(gcf,'Color',FigureBackground);
axis vis3d off;                  
lh = light('Position',[-1 1 0],'Style','infinite');
material([Ambient Diffuse Specular SpecExp SpecCol]);
view(-135, 20);




% %% =========================== PLOT SURFACES ===============================
% brainh = patch('vertices', v.surface(layer).vertices,'faces', v.surface(layer).faces,...
%     'edgecolor',ec,'BackFaceLighting',bf,...
%     'facealpha',1,'facecolor',fc,'facelighting','phong');
% 
% set(gca,'DataAspectRatio',[1 1 1]);
% set(gcf,'Color',FigureBackground);
% axis vis3d off;            
% light;
% material([Ambient Diffuse Specular SpecExp SpecCol]);
% view(90,0);

% rotate(brainh,[0 0 1],-20);

% view(35, 35);


set(0,'DefaultLineLineSmoothing','on'); 
set(0,'DefaultPatchLineSmoothing','on');



%============== Plot slices
if SliceView >= 0
%     Xlim = xlim;
%     Ylim = ylim;
%     Zlim = zlim;
    Xlim = [0 -16];
    Ylim = [-8 -20];
    Zlim = [10 -8];

    if SliceView==0         % Axial view from above, with coronal OR sagittal slices
      	Xcoord = repmat(SlicePos,[4,1]);
        Ycoord = repmat(reshape([Ylim;Ylim],[4,1]),[1,numel(SlicePos)]);
        Zcoord = repmat([Zlim, Zlim([2,1])]',[1,numel(SlicePos)]);
        view(0,90);
        
    elseif SliceView==1   	% left hemisphere sagittal view, with coronal slices
        Xcoord = repmat([Xlim, Xlim([2,1])]',[1,numel(SlicePos)]);
        Ycoord = repmat(reshape([Ylim;Ylim],[4,1]),[1,numel(SlicePos)]);
        Zcoord = repmat(SlicePos,[4,1]);
        view(-90,0);
        
    elseif SliceView==2     % Coronal, posterior
       	Xcoord = repmat([Xlim, Xlim([2,1])]',[1,numel(SlicePos)]);
        Ycoord = repmat(SlicePos,[4,1]);
        Zcoord = repmat(reshape([Zlim;Zlim],[4,1]),[1,numel(SlicePos)]);
    	view(0,0);
        
    end
    Slices = fill3(Xcoord, Ycoord, Zcoord, SliceColor,'facecolor',SliceColor,'facealpha',SliceAlpha,'edgecolor',SliceColor,'linewidth',2);
end

%============ Save image
if Save == 1
    myaa(8);
    export_fig(sprintf('CoordinateSurfaceRender_%d.png',SliceView),'-png','-transparent','-nocrop');
elseif Save > 1
    
end
