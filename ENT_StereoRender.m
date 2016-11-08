%============================ ENT_StereoRender.m ==========================
% The function renders the 3D axes with the handle specified by input 'H'
% stereoscopically for either side-by-side or red-cyan anaglyph viewing.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function ENT_StereoRender(H, Mode)

if nargin == 0
    H       = gca;
    Mode    = 1;
end
    

%=================== Set parameters
HighRes         = 0;                                        
[Az, El]        = view;                                   	% Get current azimuth and elevation angles
parallaxAngle   = [-0.5, 0.5];                             	% 
cameraAnchor    = [max(get(H, 'Xlim')), max(get(H,'Ylim')), min(get(H,'Zlim'))]; 	% Point in space to lock camera to  
set(gca, 'Units', 'pixels');
AxesPos         = get(gca, 'position');
Rect            = AxesPos;
tempfile        = 'temp.png';

%=================== Capture two views
for Eye = 1:2
%     axh(Eye) = subplot(1,2,Eye);
    axis vis3d;                                             % Lock the aspect ration for rotation
%     camtarget('auto')
    camtarget(gca, [cameraAnchor])                       	% Point the camera at the anchor point
    view(Az+parallaxAngle(Eye), El);                        % Rotate the figure to the correct angle
    if HighRes == 0                                         % For screen resolution renders...
        temp = getframe(gcf, Rect);                        	% Campture figure as image
        im{Eye} = temp.cdata;                               
    elseif HighRes > 0                                      % For high resolution redners...
        temp = export_fig(tempfile,'-png','-nocrop','-m2'); % Render each eye's image as file
        im{Eye} = imread(tempfile);                         % Read image from file
    end
end
% Link = linkprop(axh, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
% setappdata(gcf, 'StoreTheLink', Link);

if exist(tempfile, 'file')
    delete(tempfile);
end

if Mode == 1    %=================== Create anaglyph
    im{1}(:,:,2:3)  = 0 ;                                       % Removes green and blue from the left eye image
    im{2}(:,:,1)    = 0 ;                                       % Removes red from the right eye image
    stereoimage = im{1} + im{2} ;                            	% Combines the two to produce the finished anaglyph
elseif Mode == 0
    stereoimage = [im{1}, im{2}];                             	% Combines the two to produce the side-by-side image
end
figure('name','Stereo Render');                                 % Open new figure window
imshow(stereoimage,'border','tight') ;                          % Show the stereoscopic image


end