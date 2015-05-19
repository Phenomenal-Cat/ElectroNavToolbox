function [ContactCoords, SessionParams] = EN_GetContactCoordinates(Dates, SubjectDir)

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
%   SubjectDir:	the full path of the subject directory (/Subjects/SubjectID/)
%               containing: 
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
%                   - d is the number of dates requested, 
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
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%========================================================================== 


%% =========================== CHECK INPUTS AND FILES
if ~exist('SubjectDir','var') || isempty(SubjectDir)
    SubjectDir = uigetdir('','Select a subject directory');
end
RootDir = fileparts(mfilename('fullpath'));
if ~isdir(SubjectDir)
    error('No directory found for ''%s''!', SubjectDir);
end
HistoryFile = wildcardsearch(SubjectDir,'*.xls');                   
if isempty(HistoryFile)
    HistoryFile = wildcardsearch(SubjectDir,'*.csv');
end
HistoryFile = HistoryFile{1};
if isempty(HistoryFile)
    error('No recording history file (.xls/.csv) was found in ''%s''!', SubjectDir);
end
TformFile = wildcardsearch(SubjectDir,'*.mat');
if ~isempty(TformFile)
    load(TformFile{1});
end
if ~exist('T','var')  
    fprintf(['\nWARNING: \tno transformation matrix (.mat) was found in %s!\n',...
    '\t\tReturned coordinates with be in grid-centred space!\n'], SubjectDir);
end
if ~isempty(Dates)
    if numel(Dates{1})==8 && strcmp(Dates{1}(1:2),'20')             	% If input was format YYYYMMDD...
        for d = 1:numel(Dates)
            Dates{d} = YMD2DMY(Dates{d});                           	% Convert to DD-MMM-YYYY format
        end
    end
%     Dates = datestr(sort(datenum(Dates)));                           	% Check date formats and sort in chronological order
    Dates = datestr(datenum(Dates));
    SessionParams = EN_LoadSessionParams(HistoryFile, cellstr(Dates)); 	% Get electrode location information
else
    SessionParams = EN_LoadSessionParams(HistoryFile);                  % Get electrode location information
end


%% =========================== LOAD RECORDING HISTORY DATA     
GridCoords = nan(numel(SessionParams),3,24);                            % Pre-allocate at matrix to store all contact coordinates
for d = 1:numel(SessionParams)                                          % For each session date requested...
    Electrode = GetElectrodeParams(SessionParams(d).ElectrodeID);       % Get the electrode paramaters for the electrode used in this session
    if SessionParams(d).Depth > 0                                       % If tip depth (mm) is positive
        SessionParams(d).Depth = -SessionParams(d).Depth;               % Invert the depth
    end
    for c = 1:Electrode.ContactNumber                                   % For each electrode contact...
        GridCoords(d,[1,2],c) = SessionParams(d).Target;                % Get the grid hole coordinate
        GridCoords(d,3,c) = SessionParams(d).Depth+Electrode.TipLength+((c-1)*Electrode.ContactSpacing);
    end
    if exist('T','var')                                                 % If transformation matrix was loaded...
        GridCoords(d,4,:) = 1;                                          % Pad 4th row with 1s
        temp = T*squeeze(GridCoords(d,:,:));                            % Apply transformation
        ContactCoords(d,:,:) = temp(1:3,:);                             % Remove 4th row
    else
        ContactCoords(d,:,:) = GridCoords(d,:,:);
    end
end
