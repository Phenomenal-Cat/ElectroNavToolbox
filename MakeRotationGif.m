function ExportRotationMovie(fh, Filename, Duration)

%========================= ExportRotationMovie.m ==========================
% This function will create a movie file of one full rotation of a 3D plot
% about the z-axis. The scene camera angle, lighting and material settings
% should all be applied prior to calling this function. The following file
% formats are supported: .avi/.mp4/.mov/.gif.
%
% INPUTS:
%   fh:         figure handle (optional)
%   Filename:   File to save movie to, including file format extension
%   Duration:   Time in seconds for complete rotation
%   Velocity:   velocity profile, 1 = linear; 2 = sinusoidal;
%
% REFERENCES:
%   http://www.mathworks.com/help/matlab/ref/videowriter-class.html
%   http://horchler.github.io/QTWriter/
%
% 
%==========================================================================


% if exist(fh,'var')==0                               % If no figure handle was supplied...
    fh = gcf;                                       % Get handle to current figure
% end
Rect = get(fh,'position')-[0 0 0 40];               % Get the size of the figure window
% set(fh, 'toolbar','none','menu','none');            % Remove menu and toolbar


Filename = 'Neuromaps_structures.mov';

if isempty(strfind(Filename,'.'))                   % If no extension was provided...
    error('File format extension not specified in filename ''%s''!\n', Filename);
end
FileFormat = Filename(end-3:end);
switch FileFormat
    case '.avi'
        writerObj = VideoWriter(Filename, 'AVI');
    case '.mp4'
        writerObj = VideoWriter(Filename, 'MPEG-4');
    case '.mov'
        Compression = 'Photo PNG';
        Transparency = true;
        movObj = QTWriter(Filename, 'MovieFormat', Compression, 'Transparency', Transparency);
     	movObj.FrameRate = 60;  
    case '.gif'
        
    otherwise
        error('File type ''%s'' cannot be exported!\n', Filename(end-3:end));
end


Duration = 8;                                       % Duration of animation
FPS = 60;                                           % Default frame rate is 60Hz
TotalFrames = Duration*FPS;                         % Number of frames to render
                      

El = 10;
Az = 90;
Theta = linspace(Az, Az+360, TotalFrames+1);
Frame = getframe(fh, Rect);


%% ========================= RUN ANIMATION LOOP ===========================
for f = 1:TotalFrames
    view(Theta(f), El);                             % Udpate camera position
    lh = camlight('headlight');
%     lh = light('Position',[-1 1 0],'Style','infinite');
%     zoom off
%     rotate(brainh, [0 0 1], DegPerFrame);
    drawnow;                                        % Update figure
%   sfh = myaa(8);                                 	% Smooth image?
    Frame = getframe(fh, Rect);                     % Capture figure frame
                  
    switch FileFormat
        case '.mp4'     %========================== MPEG-4
            
            
        case '.mov'     %========================== QUICKTIME MOVIE                     
            writeMovie(movObj,Frame);  
            
        case '.gif'  	%========================== GIF
        	im = frame2im(Frame);                           
            [imind,cm] = rgb2ind(im,256); 
            if f==1
                imwrite(imind,cm,Filename,'gif','DelayTime',0,'loopcount',inf,'TransparentColor',1);
            elseif f <= TotalFrames
                imwrite(imind,cm,Filename,'gif','DelayTime',0,'writemode','append','TransparentColor',1);
            end
            
    end
    delete(lh);
end


switch FileFormat
    case '.avi'
        writerObj = VideoWriter(Filename, 'AVI');
    case '.mp4'
        writerObj = VideoWriter(Filename, 'MPEG-4');
    case '.mov'
        movObj.Loop = 'loop';
      	close(movObj);
    case '.gif'
        
end


