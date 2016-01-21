function MakeRotationGif(fh, Filename, Options)

%========================= ExportRotationMovie.m ==========================
% This function will create a movie file of one full rotation of a 3D plot
% about the z-axis. The scene camera angle, lighting and material settings
% should all be applied prior to calling this function. The following file
% formats are supported: .avi/.mp4/.mov/.gif.
%
% INPUTS:
%   fh:         figure handle (optional)
%   Filename:   File to save movie to, including file format extension
%   Options.Duration:   Time in seconds for complete rotation
%   Options.Velocity:   velocity profile, 1 = linear; 2 = sinusoidal;
%   Options
%
% REFERENCES:
%   http://www.mathworks.com/help/matlab/ref/videowriter-class.html
%   http://horchler.github.io/QTWriter/
%
% 
%==========================================================================

if nargin == 0
    fh                	= gcf;                           	% Get handle to current figure
    Filename            = '329_ECoG_3D.mov';
    Options.Duration    = 6;                                % Duration of animation (seconds)
    Options.Stereo      = 0;                                % 1 = red-green anaglyph; 2 = side-by-side
end

Rect = get(fh,'position');%-[0 0 0 40];                   % Get the size of the figure window
set(gca,'units','pixels');
Rect = get(gca,'position');
% set(fh, 'toolbar','none','menu','none');            % Remove menu and toolbar


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


FPS         = 60;                                               % Default frame rate is 60Hz
TotalFrames = round(Options.Duration*FPS);                     	% Number of frames to render
[Az, El]    = view;                                 
Theta       = linspace(Az, Az+360, TotalFrames+1);  
EyeOffsetDeg = -0.5;                                             % Set rotation (degrees) between two eye's perspectives
if Options.Stereo == 2
    fhsbs = figure('units','pixels','position',[1 1 1920 1080]);
%     StereoRect = get(fhsbs,'position')
    StereoRect = [1 1 1920 1080];
end

%% ========================= RUN ANIMATION LOOP ===========================
for f = 1:TotalFrames
    figure(fh);
  	view(Theta(f), El);                             % Udpate camera position
    lh = camlight('headlight');                     % Move light?
    if Options.Stereo == 0                          
        drawnow;                                        % Update figure
    %   sfh = myaa(8);                                 	% Smooth image?
        Frame = getframe(fh, Rect);                     % Capture figure frame
    elseif Options.Stereo == 2                  
        for Eye = 1:2	
            view(Theta(f)+(round(Eye-1.5)*EyeOffsetDeg), El); 	% Udpate camera position
            drawnow;
            StereoFrame{Eye} = frame2im(getframe(fh, Rect));   	% Capture figure frame
        end
        figure(fhsbs);
        set(fhsbs, 'name', sprintf('Frame %d/ %d)...', f, TotalFrames));
        FrameIm = [StereoFrame{1}, StereoFrame{2}];
        image(FrameIm);
        axis equal tight off;
        set(gca, 'units', 'pixels');
        StereoRect = get(gca, 'position');
        Frame = getframe(fhsbs, StereoRect);                  	% Capture figure frame
    end
%     lh = light('Position',[-1 1 0],'Style','infinite');
%     zoom off
%     rotate(brainh, [0 0 1], DegPerFrame);
        
    switch FileFormat
        case '.mp4'     %========================== MPEG-4
            
            
        case '.mov'     %========================== QUICKTIME MOVIE                     
            writeMovie(movObj, Frame);  
            
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


