%============================ SetDefaultParameters ========================
% This function allows the user to set the default directory and file paths
% for files required by ElectroNav. These parameters are saved to a .mat
% file in the ElectroNav directory and loaded each time the GUI is run.
%
%     ___  ______  __   __
%    /   ||  __  \|  \ |  \    APM SUBFUNCTIONS
%   / /| || |__/ /|   \|   \   Aidan P. Murphy - murphyap@mail.nih.gov
%  / __  ||  ___/ | |\   |\ \  Section of Cognitive Neurophysiology and Imaging
% /_/  |_||_|     |_| \__| \_\ Laboratory of Neuropsychology, NIMH
%==========================================================================

function Success = SetDefaultParams(Filename)

if nargin == 0
    DefaultFilename = 'ElectroNavDefaultParams.mat';
    [FileName, FilePath] = uiputfile(DefaultFilename, 'Save default parameter file as...');
    if FileName == 0
        Success = 0;
        return;
    end
    Filename = fullfile(FilePath,FileName);
end
if exist(DefaultFilename,'file')~= 0
    load(DefaultFilename);
    
end

%========================= OPEN GUI WINDOW ================================
Fig.Background = [0.5 0.5 0.5];
Fig.FontSize = 12;
Fig.TitleFontSize = 14;
Fig.PanelYdim = 130;
Fig.Rect = [400 400 500 500];                                           % Specify figure window rectangle
Fig.Handle = figure('Name','ElectroNav: Set default parameters',...   	% Open a figure window with specified title
                    'Color',Fig.Background,...                          % Set the figure window background color
                    'Renderer','OpenGL',...                             % Use OpenGL renderer
                    'OuterPosition', Fig.Rect,...                       % position figure window
                    'NumberTitle','off',...                             % Remove figure number from title
                    'Menu','none','Toolbar','none');                    % Turn off toolbars to save space
Fig.Margin = 20;
BoxPos{1} = [Fig.Margin, Fig.Rect(4)-160, Fig.Rect(3)-(2*Fig.Margin), Fig.PanelYdim];   
BoxPos{2} = [Fig.Margin, BoxPos{1}(2)-BoxPos{1}(4)-Fig.Margin, Fig.Rect(3)-(2*Fig.Margin), Fig.PanelYdim]; 
BoxPos{3} = [Fig.Margin, BoxPos{2}(2)-BoxPos{2}(4)-Fig.Margin, Fig.Rect(3)-(2*Fig.Margin), Fig.PanelYdim];
LabelDim = [420, 20];
BoxDim = [80 20];
ButtonOptions = {'Select','Skip','Info'};
ButtonCallbacks = {'GetFile', 'SkipFile','FileInfo'};

%======================= Physiology data
PhysioLabels = {'Select Microsoft Excel file (*.xls; *.xlsx) containing record of physiology recording sessions',...
                'Select directory containing analysed physiology data for all subjects'};
PhysioHandle = uipanel( 'Title','Electrophysiology Data',...
                        'FontSize',Fig.TitleFontSize,...
                        'BackgroundColor',Fig.Background,...
                        'Units','pixels',...
                        'Position',BoxPos{1},...
                        'Parent',Fig.Handle);

for i = 1:numel(PhysioLabels)
    Pos = numel(PhysioLabels)-i;
    PhysioLabelPos{i} = [10, 10+Pos*LabelDim(2),LabelDim];
    PhysioLabelHandle(i) = uicontrol( 'Style','Text',...
                                        'String',PhysioLabels{i},...
                                        'FontSize',Fig.FontSize,...
                                        'HorizontalAlignment','Left',...
                                        'pos',PhysioLabelPos{i},...
                                        'BackgroundColor',Fig.Background,...
                                        'parent',PhysioHandle);
 	for b = 1:numel(ButtonOptions)
        PhysioInputPos = [10+(b-1)*(BoxDim(1)+10),PhysioLabelPos{i}(2)-BoxDim(2),BoxDim];
        PhysioInputHandle(i+b) = uicontrol( 'Style','Pushbutton',...
                                            'String',ButtonOptions{b},...
                                            'HorizontalAlignment','Left',...
                                            'pos',PhysioInputPos,...
                                            'parent',PhysioHandle);
%                                             'Callback',{@ButtonCallbacks{b},i},...

    end
end

%======================= MRI data
MRILabels = {   'Select Nifti file (*.nii; *.img & *.hdr) containing anatomical grid scan',...
                'Select Nifti file (*.nii; *.img & *.hdr) containing atlas volume in native space',...
                'Select directory containing MRI data data for all subjects'};
MRIDetails = {'Select','Skip'};
MRIHandle = uipanel('Title','MRI Data',...
                    'FontSize',Fig.TitleFontSize,...
                    'BackgroundColor',Fig.Background,...
                    'Units','pixels',...
                    'Position',BoxPos{2},...
                    'Parent',Fig.Handle);
                
for i = 1:numel(MRILabels)
    Pos = numel(MRILabels)-i;
    MRILabelPos{i} = [10, 10+Pos*LabelDim(2),LabelDim];
    MRILabelHandle(i) = uicontrol( 'Style','Text',...
                                        'String',MRILabels{i},...
                                        'HorizontalAlignment','Left',...
                                        'FontSize',Fig.FontSize,...
                                        'pos',MRILabelPos{i},...
                                        'BackgroundColor',Fig.Background,...
                                        'parent',MRIHandle);
    MRIInputHandle(i) = uicontrol( 'Style','Pushbutton',...
                                        'String',MRIDetails{i},...
                                        'HorizontalAlignment','Left',...
                                        'pos',[BoxDim(1)+20,10+Pos*BoxDim(2),80,20],...
                                        'parent',MRIHandle);
end


%======================= 3D surface data
SurfLabels = {'Select X-form matrix file (*.xform; *.txt) containing grid -> MRI X-form transformation matrix',...
                'Select directory containing 3D surface data (*.vtk) for all structures'};
SurfDetails = {'Select','Skip'};
SurfHandle = uipanel(   'Title','Surf Data',...
                        'FontSize',Fig.TitleFontSize,...
                        'BackgroundColor',Fig.Background,...
                        'Units','pixels',...
                        'Position',BoxPos{3},...
                        'Parent',Fig.Handle);
                    
for i = 1:numel(SurfLabels)
    Pos = numel(SurfLabels)-i;
    SurfLabelPos{i} = [10, 10+Pos*LabelDim(2),LabelDim];
    SurfLabelHandle(i) = uicontrol( 'Style','Text',...
                                        'String',SurfLabels{i},...
                                        'HorizontalAlignment','Left',...
                                        'FontSize',Fig.FontSize,...
                                        'pos',SurfLabelPos{i},...
                                        'BackgroundColor',Fig.Background,...
                                        'parent',SurfHandle);
    SurfInputHandle(i) = uicontrol( 'Style','Pushbutton',...
                                        'String',SurfDetails{i},...
                                        'HorizontalAlignment','Left',...
                                        'pos',[BoxDim(1)+20,10+Pos*BoxDim(2),80,20],...
                                        'parent',SurfHandle);
end










% Defaults.AllDirectories = {	Defaults.Files.ExcelFile,  ...
%                         Defaults.Files.PhysioRootDir, ...
%                         Defaults.Files.MRIRootDir,  ...
%                         Defaults.Files.MRImage,  ...
%                         Defaults.Files.Xformfile,...
%                         Defaults.Files.VTKfile};


    
    
    
    
    
    
    
    
    
end


function GetFileh(Obj, Event, Indx)
    disp('GetFile')


end

function SkipFile(hObj, Event, Indx)
    disp('SkipFile')

end

