%============================= EN_About.m =================================
% This subfunction presents the ElectroNav logo and toolbox information in
% a figure window, either during loading or on request. It returns the
% figure handle as an output so that it can be closed once loading has
% completed.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy (murphyap@mail.nih.gov)
% © Copyleft 2014, GNU General Public License
%==========================================================================

function FigHandle = EN_About(Loading, temp)
    if nargin ==0
        Loading = 0;
    end
    ElectroNavVersion = 1.0;
  	ContactEmail = 'murphyap@mail.nih.gov';
    GithubSite = 'https://github.com/MonkeyGone2Heaven/ElectroNavToolbox';
    FigBackground = [1 1 1]*0.75;
    LogoFile = 'ElectroNavLogo1.png';
    FigLogo = imread(LogoFile,'BackgroundColor',FigBackground);
    [a b FigLogoAlpha] = imread(LogoFile); 
    LogoSize = size(FigLogo);
    FigAboutRect = [500 500 500 300];                
    FigMargin = 40;
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
    LogoAxH = axes('Units','pixels','position',[FigMargin,160,LogoSize(2),LogoSize(1)],'visible','off');
    LogoH = image(FigLogo);                                         % Display the logo image
    set(LogoH, 'AlphaData',FigLogoAlpha)
    axis equal tight off;                                           % turn axes off
    
    %============= Add text and links to figure window
    Text1 = 'A MATLAB® Toolbox for MRI-guided Electrode Navigation';
    Text2 = sprintf(['Version %.1f, developed by Aidan Murphy %s Copyleft 2014\n',...
                    'Section on Cognitive Neurophysiology and Imaging, NIMH\n\n',...
                    'Contact:\n\nDownload:'],ElectroNavVersion, char(169));
  	TextAxH(1) = axes('Units','pixels','position',[FigMargin,FigMargin-10,FigAboutRect(3)-(2*FigMargin),FigAboutRect(4)*0.38],'visible','off');
  	TextH(1) = text(0,1,Text1,'FontWeight','bold','FontUnits','pixels','FontSize',15,'HorizontalAlignment','left','VerticalAlignment','top');
    TextAxH(2) = axes('Units','pixels','position',[FigMargin,40,FigAboutRect(3)-(2*FigMargin),FigAboutRect(4)*0.18],'visible','off');
    TextH(2) = text(0,1,Text2,'FontUnits','points','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top');
    labelStr{1} = sprintf('<html><center><a href="">%s', ContactEmail);                     % Set label for e-mail button
    labelStr{2} = sprintf('<html><center><a href="">GitHub'); 
    cbStr{1} = sprintf('web([''mailto:'',''%s'']);', ContactEmail);                       	% Set link to e-mail
    cbStr{2} = sprintf('web(''%s'');', GithubSite);                                             
    hButton(1) = uicontrol('string',labelStr{1},'pos',[110,45,160,20],'callback',cbStr{1});	% Create push button to send e-mail
    hButton(2) = uicontrol('string',labelStr{2},'pos',[110,20,160,20],'callback',cbStr{2}); 
    
    if Loading == 1                                                                         % If input 'Loading' is 1...
        TextAxH(3) = axes('Units','pixels','position',[FigMargin,FigAboutRect(4)-(2*FigMargin),500,40],'visible','off');
        TextH(3) = text(0,1,'Loading...','FontWeight','bold','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top');
    end
    drawnow;                                                                                % Display figure window
end