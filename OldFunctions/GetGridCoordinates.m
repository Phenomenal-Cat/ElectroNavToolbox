function [GridCoords, Dates] = GetGridCoordinates(Dates, Excel, Sheet)

%========================= GetGridCoordinates.m ===========================
% This file reads the specified Excel file, containing electrode positions
% for each past recording session, and returns the 3D coordinates of all
% electrode contact positions in a grid-centred reference frame (the centre 
% grid hole at the base of the grid is the origin).
%
% INPUTS:
%   Dates:      a 1 x d cell array containing date strings in the format 
%               'DD-MMM-YYYY' or 'YYYYMMDD'. Leaving Dates input empty
%               allows the user to manually select dates from a list of all
%               dates in the Excel file.
%   Excel:      an Excel file containing electrode position data for all
%               session dates.
%   Sheet:      which number spreadsheet to read data from.
%
% OUTPUTS:
%   GridCoords: a d x c x n matrix containing grid-centred coordinates in
%               mm, where d is the number of dates requested, c is x, y and 
%               z coordinates and n is the number of electrode contacts for 
%               the electrode used. Empty cells are filled with NaNs.
%   Dates:      optional. If the input Dates was empty, a cell array of
%               manually selected dates is returned.
%
% REVISIONS:
%   12/10/2014 - Written by Aidan Murphy
%   16/10/2014 - Manual date selection added
%     ___  ______  __   __
%    /   ||  __  \|  \ |  \    APM SUBFUNCTIONS
%   / /| || |__/ /|   \|   \   Aidan P. Murphy - murphyap@mail.nih.gov
%  / __  ||  ___/ | |\   |\ \  Section of Cognitive Neurophysiology and Imaging
% /_/  |_||_|     |_| \__| \_\ Laboratory of Neuropsychology, NIMH
%==========================================================================


%=========================== CHECK INPUTS
if ~exist('Excel','var')        
    error('Excel file input was not provided!\n');
end
if ~exist(Excel,'file')
    error('Excel file ''%s'' does not exist!\n', Excel);
end
if ~strcmpi(Excel(end-2:end),'xls')                                 	% If selected history file was not Excel format...
    error('Excel input ''%s'' is not .xls format!\n', Excel);
end
if ~isempty(Dates)
    if numel(Dates{1})==8                                            	% If input was format YYYYMMDD...
        for d = 1:numel(Dates)
            Dates{d} = YMD2DMY(Dates{d});                           	% Convert to DD-MMM-YYYY format
        end
    end
%     Dates = datestr(sort(datenum(Dates)));                           	% Check date formats and sort in chronological order
    Dates = datestr(datenum(Dates));
end

%=========================== LOAD EXCEL DATA
[status,SheetNames] = xlsfinfo(Excel);                                	% Get Excel sheet names
if ~exist('Sheet','var')
    [Sheet,ok] = listdlg('ListString',SheetNames,...                	% Ask user to select a sheet
                         'ListSize',[160 60],...
                         'SelectionMode', 'multiple',...
                         'PromptString','Select Excel sheet(s):'); 
    if ok==0, return; end
end
Data = [];
for s = 1:numel(Sheet)
    [num,txt,raw] =  xlsread(Excel, SheetNames{Sheet(s)},'','basic'); 	% Read Excel file
    Data = [Data; num];
end
Headers = txt{1,:};
if isnan(num(1,1))
    num(1,:) = [];                                                     	% Remove nans
end
DateNums = num(:,1)+datenum('30-Dec-1899');                            	% Convert Excel dates to Matlab dates
DateStrings = datestr(DateNums);                                                                  

%========================= MANUALLY SELECT FROM AVAILABLE DATES
if isempty(Dates)
     [Selection,ok] = listdlg('ListString',DateStrings,'SelectionMode','multi','PromptString','Select date(s) of previous session:');
    if ok==0, return; end
    Dates = cellstr(datestr(DateNums(Selection)));
end

%========================== GET ELECTRODE TIP COORDINATES FOR SPECIFIED DATES
for d = 1:numel(Dates(:,1))
    SessionIndx = find(DateNums == datenum(Dates(d,:)));
    if isempty(SessionIndx)
        error('Input date ''%s'' was not found in %s!\n', Dates(d,:), Excel);
    else
        SessionIndices(d) = SessionIndx;
    end
end
TipCoords = num(SessionIndices,[2,3,4]);
TipCoords(:,3) = -TipCoords(:,3);                                     % Convert depth coordinate to negative (mm relative to grid base)

%==================== CACLULATE POSITION OF EACH CONTACT
GridCoords = nan(numel(SessionIndices),3,24);
ElectrodeType = txt(SessionIndices+1,7);
for e = 1:numel(ElectrodeType)
    Electrode = GetElectrodeParams(ElectrodeType{e});
    for c = 1:Electrode.ContactNumber
        GridCoords(e,[1,2],c) = TipCoords(e,[1,2]);
        GridCoords(e,3,c) = TipCoords(e,3)+Electrode.TipLength+((c-1)*Electrode.ContactSpacing);
    end
end
