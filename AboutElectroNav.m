%============================ AboutElectroNav =============================
% This subfunction presents the ElectroNav logo and toolbox information in
% a figure window, either during loading or on request. It returns the
% figure handle as an output so that it can be closed once loading has
% completed.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy (murphyap@mail.nih.gov)
% © Copyleft 2014, GNU General Public License
%==========================================================================

function FigHandle = AboutElectroNav(Loading)
    if nargin ==0
        Loading = 0;
    end
 	BackgroundOn = 0;
    ElectroNavVersion = 1.0;
    FigBackground = [1 1 1]*0.75;
    LogoFile = 'ElectroNav_5.png';
    FigLogo = imread(LogoFile,'BackgroundColor',FigBackground);
    [a b FigLogoAlpha] = imread(LogoFile); 
    LogoSize = size(FigLogo);
    ContactEmail = 'murphyap@mail.nih.gov';
    FigAboutRect = [500 500 500 300];                
    FigMargin = 40;
%    	fdr=@(varargin) 'disp(''DELETE'');';
%      fcr=@(varargin) 'disp(''CLOSE'');';
%      fcf=@(varargin) set(gcbo,...                                   % Set callback for figure creation
%                     'closerequestfcn',fcr(),...                     
%                     'deletefcn',fdr(),...                           
%                     'handlevisibility','off');                      %
    fcf = [];
    FigHandle = figure('Name',['ElectroNav',char(169)],...        	% Open a figure window with specified title
                        'Color',FigBackground,...                   % Set the figure window background color
                        'Renderer','OpenGL',...                     % Use OpenGL renderer
                        'OuterPosition', FigAboutRect,...           % position figure window
                        'NumberTitle','off',...                     % Remove the figure number from the title
                        'Resize','off',...                          % Turn off resizing of figure
                        'createfcn',fcf,...                         % Set create callback
                        'Menu','none','Toolbar','none');            % Turn off toolbars and menu
                    
   	%============= Display logo image
    if BackgroundOn == 1
        BackgFile = 'StereoCortex.png';
        FigBackg = imread(BackgFile);
        BackgAxH = axes('Units','pixels','position',[0,0,FigAboutRect([3,4])],'visible','off');     
        image(FigBackg);                
        axis equal tight off;                                    	% turn axes off
    end
    LogoAxH = axes('Units','pixels','position',[FigMargin,160,LogoSize(2),LogoSize(1)],'visible','off');
    LogoH = image(FigLogo);                                         % Display the logo image
    set(LogoH, 'AlphaData',FigLogoAlpha)
    axis equal tight off;                                           % turn axes off
    
    %============= Add text to figure window
    Text1 = 'A MATLAB® Toolbox for MRI-guided Electrode Navigation';
    Text2 = sprintf(['Version %.1f, developed by Aidan Murphy %s Copyleft 2014\n',...
                    'Section on Cognitive Neurophysiology and Imaging, NIMH\n',...
                    'Contact:'],ElectroNavVersion, char(169));
  	TextAxH(1) = axes('Units','pixels','position',[FigMargin,FigMargin,FigAboutRect(3)-(2*FigMargin),FigAboutRect(4)*0.38],'visible','off');
  	TextH(1) = text(0,1,Text1,'FontWeight','bold','FontUnits','pixels','FontSize',15,'HorizontalAlignment','left','VerticalAlignment','top');
    TextAxH(2) = axes('Units','pixels','position',[FigMargin,20,FigAboutRect(3)-(2*FigMargin),FigAboutRect(4)*0.18],'visible','off');
    TextH(2) = text(0,1,Text2,'FontUnits','points','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top');

%     labelStr = sprintf('<html><center><a href="">%s',ContactEmail);
    labelStr = sprintf('<html><center><a href="">murphyap@mail.nih.gov');               % Set label for e-mail button
    cbStr =  'web([''mailto:'',''murphyap@mail.nih.gov'']);';                           % Set link to e-mail
    hButton = uicontrol('string',labelStr,'pos',[90,25,160,20],'callback',cbStr);       % Create push button to send e-mail
    
    if Loading == 1                                                                     % If input 'Loading' is 1...
        TextAxH(3) = axes('Units','pixels','position',[FigMargin,FigAboutRect(4)-(2*FigMargin),500,40],'visible','off');
        TextH(3) = text(0,1,'Loading...','FontWeight','bold','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top');
    end
end