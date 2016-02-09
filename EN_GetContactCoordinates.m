function [ContactCoords, SessionParams] = EN_GetContactCoordinates(Dates, SubjectID)

%======================== EN_GetContactCoordinates.m ======================
% This file reads the recording history file for the specified subject, which
% contains grid coordinates and electrode depth for each past recording 
% session. The function returns the 3D coordinates of all electrode contact 
% positions in: 
%
%       1)  a world-centred reference frame, where the anterior comissure or 
%           the interaural line is the origin (depending on the alignment of
%           the MRI volume from which the transform matrix was derived).
%       2)  If no transfromation matrix file is found then a warning appears 
%           and grid-centred coordinates are returned, where the centre grid 
%           hole at the base of the grid is the origin.
%
% INPUTS:
%   Dates:      a 1 x d cell array containing date strings in the format 
%               'DD-MMM-YYYY' or 'YYYYMMDD'. Leaving Dates input empty
%               allows the user to manually select dates from a list of all
%               dates in the Excel file.
%   SubjectID:	string containing subject ID. This will be used to load the
%               default parameters saved for the specified subject, including:
%               1)  an Excel or .cvs file containing all electrode position 
%                   data for all session dates.
%               2)  a .mat file containing the 4x4 transformation matrix for
%                   converting grid-centred coordinates into real-world
%                   coordinates (mm relative to the anterior comissure in an
%                   ACPC aligned volume).
%
% OUTPUTS:
%   ContactCoords: [d, 3, n] matrix containing world-centred coordinates in
%               mm relative to the anterior comissure, where:
%                   - d is the number of dates/ penetrations requested, 
%                   - x, y and z coordinates 
%                   - n is the number of electrode contacts for the electrode 
%               used. Empty cells are filled with NaNs.
%   SessionParams:  A structure containing parameters for each session date
%               	requested.
%
% REVISIONS:
%   16/10/2014 - Manual date selection added
%   17/03/2015 - Simplified and grid -> ACPC coordinate transform added
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, � Copyleft 2014-2016, GNU General Public License
%========================================================================== 


%% =========================== CHECK INPUTS AND LOCATE FILES
if ~exist('SubjectID','var')
    SubjectID = [];
end

Defaults    = ENT_LoadDefaults(SubjectID);
if exist(Defaults.Xform, 'file')
    load(Defaults.Xform);
else
    error('Default transform matrix file ''%s'' does not exist!', Defaults.Xform)
end

if ~exist('T','var')  
    fprintf(['\nWARNING: \tno transformation matrix (.mat) was found in %s!\n',...
    '\t\tReturned coordinates with be in grid-centred space!\n'], SubjectDir);
end
if exist('Dates','var') && ~isempty(Dates)
    if ~iscell(Dates)                                                   % If dates were not provided in a cell array...
        Dates = cellstr(Dates);                                       	% Convert to cell array
    end
    if numel(Dates{1})==8 && strcmp(Dates{1}(1:2),'20')             	% If input was format YYYYMMDD...
        for d = 1:numel(Dates)
            Dates{d} = YMD2DMY(Dates{d});                           	% Convert to DD-MMM-YYYY format
        end
    end
%     Dates = datestr(sort(datenum(Dates)));                           	% Check date formats and sort in chronological order
    Dates = datestr(datenum(Dates));
    SessionParams = ENT_LoadSessionParams(Defaults.HistoryFile, cellstr(Dates)); ...
    % Get electrode location information
else
    SessionParams = ENT_LoadSessionParams(Defaults.HistoryFile);               	% Get electrode location information
end


%% =========================== LOAD RECORDING HISTORY DATA     
NoContacts  = 24;                                                               
GridCoords  = nan(numel(SessionParams),3,NoContacts*2);                       	% Pre-allocate at matrix to store all contact coordinates
for d = 1:numel(SessionParams)                                                  % For each session date requested...
    Electrode = cell2struct(SessionParams(d).ElectrodeID,'ID');                 % Convert electrode ID(s) to structure
    Electrode = ENT_GetElectrodeParams(Electrode);                              % Get the electrode paramaters for the electrode(s) used in this session
    ChIndx = 1;                                                                 % Keep count of channel number (includes all electrodes within a session)
    for e = 1:numel(Electrode)                                                  % For each electrode used...
        if SessionParams(d).Depth{e} > 0                                        % If tip depth (mm) is positive
            SessionParams(d).Depth{e} = -SessionParams(d).Depth{e};          	% Invert the depth
        end
        for c = 1:Electrode(e).ContactNumber                                    % For each electrode contact...
            GridCoords(d,[1,2],ChIndx) = SessionParams(d).Target{e};          	% Get the grid hole coordinates
            GridCoords(d,3,ChIndx) = SessionParams(d).Depth{e}+Electrode(e).TipLength+((c-1)*Electrode(e).ContactSpacing);
         	ChIndx = ChIndx+1;
        end
    end
    
    if exist('T','var')                                                     % If transformation matrix was loaded...
        GridCoords(d,4,:) = 1;                                              % Pad 4th row with 1s
        temp = T*squeeze(GridCoords(d,:,:));                                % Apply transformation
        ContactCoords(d,:,1:size(temp,2)) = temp(1:3,:);                  	% Remove 4th row
    else
        ContactCoords(d,:,:) = GridCoords(d,:,:);
    end
end



