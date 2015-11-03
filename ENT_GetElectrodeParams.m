%=========================== ENT_GetElectrodeParams =======================
% This function returns a structure ('Electrode') containing information
% about the physical properties of the electrode specified by the input
% structure field 'Electrode.ID'. These parameters are hard coded in this 
% function, based on currently available linear multielectrode array 
% technology from four major manufacturers. To add additional manufacturers 
% or to change the defaults for an existing manufacturer, edit the parameters 
% listed.
%
% ElectrodeID:  Manufacturer:                   
%               'AO'    = Alpha Omega           
%               'NN'    = NeuroNexus             
%               'PLX'   = Plexon          
%               'TR'    = Thomas Recording 
%   
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function Electrode = ENT_GetElectrodeParams(Electrode)

%================ CHECK INPUTS
Brands = {'AO','NN','PLX','TR'};
LineNumber = 22;                                                            % Line number to edit from
if nargin == 0                                                              % If no input variable...
    Electrode = Brands;                                                   	% Return list of accepted inputs
    return;                                                                 % Exit function
end
if ~isstruct(Electrode)
    ElectrodeType = Electrode;
    clear Electrode;
    Electrode.ID = ElectrodeType; 
end
NoElectrodes = length(Electrode);                                           % If Electrode structure contains data for multiple electrodes...

%================ FOR EACH ELECTRODE...
for e = 1:NoElectrodes                                                      % For each electrode...
    Electrode(e).AllTypes = Brands;
    if ~isfield(Electrode,'ID') || isempty(Electrode(e).ID)
        if isfield(Electrode,'Brand') && ~isempty(Electrode(e).Brand)
            Electrode(e).ID = Electrode(e).Brand;
        elseif ~isfield(Electrode,'Brand') || isempty(Electrode(e).Brand)
            error('Input structure ''Electrode'' must contain ''ID'' field!');
        end
    end 
    
    if ~isfield(Electrode, 'ContactNumber') || isempty(Electrode(e).ContactNumber)
        Electrode(e).ContactNumber = str2num(Electrode(e).ID(regexp(Electrode(e).ID,'\d')));                    % Extract number of contacts from ID
    end
    Electrode(e).Brand = Electrode(e).ID(1:findstr(Electrode(e).ID, num2str(Electrode(e).ContactNumber))-1);    % Extract manufacturer abbreviation from ID

    switch Electrode(e).ID(1:2)

        case 'PL'                      %============ Plexon V-probe/U-probe linear multi-electrode array
            Electrode(e).Length = 90;                  % Full electrode shaft length (mm)
            Electrode(e).Diameter = 0.2;               % shadt diameter (mm)
            Electrode(e).TipLength = 1.5;              % distance from tip to first contact (mm)
            Electrode(e).ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
            Electrode(e).ContactDiameter = 0.1;        % Exagerate contact diameter for visualization

        case 'AO'                       %============ Alpha Omega linear multi-electrode array (microfil)
            Electrode(e).Length = 70;                  % Full electrode shaft length (mm)
            Electrode(e).Diameter = 0.3;               % shadt diameter (mm)
            Electrode(e).TipLength = 1.5;              % distance from tip to first contact (mm)
            Electrode(e).ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
            Electrode(e).ContactDiameter = 0.1;        % Exagerate contact diameter for visualization   

        case 'NN'                       %============ NeuroNexus linear multi-electrode array
            Electrode(e).Length = 70;                  % Full electrode shaft length (mm)
            Electrode(e).Diameter = 0.3;               % shadt diameter (mm)
            Electrode(e).TipLength = 1;                % distance from tip to first contact (mm)
            Electrode(e).ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
            Electrode(e).ContactDiameter = 0.1;        % Exagerate contact diameter for visualization

        case 'TR'                       %============ Thomas Recording
            Electrode(e).Length = 70;                  % Full electrode shaft length (mm)
            Electrode(e).Diameter = 0.3;               % shadt diameter (mm)
            Electrode(e).TipLength = 1;                % distance from tip to first contact (mm)
            Electrode(e).ContactSpacing = 0.3;         % distance between adjacent contacts (mm)
            Electrode(e).ContactDiameter = 0.1;        % Exagerate contact diameter for visualization

        otherwise                       %============ Unknown manufacturer: suggest user edits m-file                 
            Msg = sprintf([ 'The electrode ID ''%s'' does not belong to a recognized manufacturer!\n', ...
                            'Please edit <a href="matlab: opentoline(which('' %s.m''),%d)">%s.m</a> to add parameters for a new electrode brand.'], Electrode(e).ID, mfilename, LineNumber, mfilename);
            disp(Msg);
            h = warndlg(Msg,'Electrode Error');
            return;
    end

        %=================== Default electrode visualization settings =============
    if ~isfield(Electrode, 'CurrentSelected') || isempty(Electrode(e).CurrentSelected)
        Electrode(e).CurrentSelected = 1;                                       % Default contact selection is #1 (nearest to tip)
    end
    if ~isfield(Electrode,'ContactData') || isempty(Electrode(e).ContactData)
        Electrode(e).ContactData = zeros(Electrode(e).ContactNumber,1);         % Default contact availability is none
    end
    Electrode(e).Color              = [0.4 0.4 0.4];                         	% Schematic electrode shaft color
    Electrode(e).ContactColor       = [1 0 0];                                  % Schematic electrode contact color
    Electrode(e).ContactLength      = Electrode(e).ContactSpacing*Electrode(e).ContactNumber; 
    Electrode(e).SelectionColor     = [1 1 1];                                  % Set highlight color for selected contacts
    Electrode(e).GuideColor         = [0,0,1];                                  % Set default guide tube color
    Electrode(e).GuideAlpha         = 0.5;                                      % Set default guide tube opacity
    Electrode(e).MRIColor           = [0 1 1];                                  % MRI view electrode shaft color
    Electrode(e).QualityColorMap	= [0 0 0; 1 0 0; 1 0.5 0; 1 1 0; 0 1 0];
end