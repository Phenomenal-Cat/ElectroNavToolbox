%============================ GetElectrodeParams ==========================
% This function returns a structure ('Electrode') containing information
% about the physical properties of the electrode specified by the input
% structure field 'Electrode.ID'. These parameters are hard coded in this 
% function, based on currently available linear multielectrode array 
% technology from four major manufacturers. To add additional manufacturers 
% or to change the defaults for an existing manufacturer, edit the parameters 
% listed.
%
% ElectrodeID:
%   Manufacturer:                   Contact number:
%       'AO' = Alpha Omega                  12
%       'NN' = NeuroNexus                   24
%       'PLX' = Plexon                      32
%       'TR' = Thomas Recording             64
%   
% 
%==========================================================================

function Electrode = GetElectrodeParams(Electrode)

Brands = {'PLX','AO','NN','TR'};
LineNumber = 22;                                                            % Line number to edit from
if nargin == 0                                                              % If no input variable...
    Electrode = Brands;                                                   	% Return list of accepted inputs
    return;                                                                 % Exit function
end
if ~isstruct(Electrode)
    ElectrodeType = Electrode;
    clear Electrode;
    Electrode.ID = ElectrodeType; 
else
    if ~isfield(Electrode,'ID') && isfield(Electrode,'Brand')
        Electrode.ID = Electrode.Brand;
    elseif ~isfield(Electrode,'ID') && ~isfield(Electrode,'Brand')
        error('Input structure ''Electrode'' must contain ''ID'' field!');
    end
end
if ~isfield(Electrode, 'ContactNumber')
    Electrode.ContactNumber = str2num(Electrode.ID(regexp(Electrode.ID,'\d')));                 % Extract number of contacts from ID
end
Electrode.Brand = Electrode.ID(1:findstr(Electrode.ID, num2str(Electrode.ContactNumber))-1);    % Extract manufacturer abbreviation from ID

switch Electrode.Brand
    
    case 'PLX'                      %============ Plexon V-probe/U-probe linear multi-electrode array
     	Electrode.Length = 90;                  % Full electrode shaft length (mm)
        Electrode.Diameter = 0.2;               % shadt diameter (mm)
        Electrode.TipLength = 1.5;              % distance from tip to first contact (mm)
        Electrode.ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization

    case 'AO'                       %============ Alpha Omega linear multi-electrode array (microfil)
        Electrode.Length = 70;                  % Full electrode shaft length (mm)
        Electrode.Diameter = 0.3;               % shadt diameter (mm)
        Electrode.TipLength = 1.5;              % distance from tip to first contact (mm)
        Electrode.ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization   
        
    case 'NN'                       %============ NeuroNexus linear multi-electrode array
        Electrode.Length = 70;                  % Full electrode shaft length (mm)
        Electrode.Diameter = 0.3;               % shadt diameter (mm)
        Electrode.TipLength = 1;                % distance from tip to first contact (mm)
        Electrode.ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization
        
    case 'TR'                       %============ Thomas Recording
        Electrode.Length = 70;                  % Full electrode shaft length (mm)
        Electrode.Diameter = 0.3;               % shadt diameter (mm)
        Electrode.TipLength = 1;                % distance from tip to first contact (mm)
        Electrode.ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
        Electrode.ContactDiameter = 0.1;        % Exagerate contact diameter for visualization
        
    otherwise                       %============ Unknown manufacturer: suggest user edits m-file                 
        Msg = sprintf([ 'The electrode ID ''%s'' does not belong to a recognized manufacturer!\n', ...
                        'Please edit <a href="matlab: opentoline(which('' %s.m''),%d)">%s.m</a> to add parameters for a new electrode brand.'], Electrode.ID, mfilename, LineNumber, mfilename);
        disp(Msg);
        h = warndlg(Msg,'Electrode Error');
        return;
end

%=================== Default electrode visualization settings =============
Electrode.Color = [0.4 0.4 0.4];                         	% Schematic electrode shaft color
Electrode.ContactColor = [1 0 0];                           % Schematic electrode contact color
if ~isfield(Electrode, 'CurrentSelected')
    Electrode.CurrentSelected = 1;                        	% Default contact selection is #1 (nearest to tip)
end
if ~isfield(Electrode,'ContactData')
    Electrode.ContactData = zeros(1,Electrode.ContactNumber);  % Default contact availability is none
end
Electrode.ContactLength = Electrode.ContactSpacing*Electrode.ContactNumber; 
Electrode.SelectionColor = [1 1 1];                         % Set highlight color for selected contacts
Electrode.MicrodriveDepth = 0;                              % Initial microdrive depth = 0mm
Electrode.StartDepth = 0;                                   % Set default electrode start depth (mm below grid base)
Electrode.GuideColor = [0,0,1];                             % Set default guide tube color
Electrode.GuideAlpha = 0.5;                                 % Set default guide tube opacity
Electrode.MRIColor = [0 1 1];                               % MRI view electrode shaft color