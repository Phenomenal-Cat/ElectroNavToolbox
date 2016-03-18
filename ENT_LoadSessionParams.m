function [SessionParams] = ENT_LoadSessionParams(HistoryFile, SessionDate)

%========================== ENT_LoadSessionParams.m ========================
% This function retreives the session parameters stored in the provided
% spreadsheet file, for the date specified.
%
% INPUTS:
%       HistoryFile:    full path of spreadsheet (.xls/ .csv) containing recording history
%       SessionDate:    optional string specifying which date to query
%                       parameters for, in DD-MMM-YYYY format. If not provided,
%                       user must select from a list. 
%                       - 'All' returns all dates available in file
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy
% ? Copyleft 2014, GNU General Public License
%==========================================================================

%========================== Check inputs
if nargin == 0
    [file, path] = uigetfile({'*.xls;*.csv;*.mat'},'Select recording history');
    HistoryFile = fullfile(path, file);
end
if ~exist(HistoryFile,'file')
    error('Session history file %s does not exist!', HistoryFile);
end

%========================== Load recording history data
[a,b,HistoryFormat] = fileparts(HistoryFile);
if exist('datetime.m','file')==2                   	%============ MATLAB R2014a and later
    switch HistoryFormat
        case '.xls'
            T           = readtable(HistoryFile,'basic',1);
            T.Date      = datetime(T.Date,'ConvertFrom','excel');
            DateStrings = datestr(T.Date);
            C           = table2cell(T);
            Header      = T.Properties.VariableNames;
        
        case '.csv'
    %         formatSpec = '%{dd-MMM-yyyy}D%f%f%f%f%f%C';
    %         T = readtable(HistoryFile,'Delimiter',',','Format',formatSpec);
            fid         = fopen(HistoryFile,'rt');
            Header      = textscan(fid, '%s%s%s%s%s%s%s\n', 1, 'delimiter', ',');
            data        = textscan(fid, '%f %f %f %f %f %f %s','headerlines',1,'delimiter',',');
            fclose(fid);
            Dates       = data{1};   
            DateStrings = datestr(Dates);
            C           = num2cell(cell2mat(data(1:6)));
            C(:,7)      = data{7};
            
        case '.mat'
            load(HistoryFile);
            Header      = raw(1,:);
            Dates       = raw(2:end,1)+datenum('30-Dec-1899');          % Convert Excel dates to Matlab dates
            DateStrings = datestr(Dates);      
            C           = raw(2:end,:);
            
        otherwise
            error('File ''%s'' is not a valid spreadsheet format!', HistoryFile);
    end
else                                                %============ MATLAB R2013b and earlier 
    [status, sheets]	= xlsfinfo(HistoryFile);
    [num,txt,raw]       =  xlsread(HistoryFile,sheets{1});         	% Read data from Excel file
    Header              = txt(1,:);                               	% Skip row containing column titles
    num(1,:)            = [];                                      	% Remove header NaNs
    Dates               = num(:,1)+datenum('30-Dec-1899');        	% Convert Excel dates to Matlab dates
    DateStrings         = datestr(Dates);                               
    C                   = raw(2:end,:);                                           
end


%========================== Get a date input
if ~exist('SessionDate','var')
    [Selection,ok] = listdlg('ListString',DateStrings,'SelectionMode','multiple','PromptString','Select session:');
    if ok==0
        SessionParams = [];
        return;
    end
    SessionDate = mat2cell(DateStrings(Selection,:),ones(1,numel(Selection)),size(DateStrings,2));
end

if ischar(SessionDate) && strcmpi(SessionDate, 'all')
    SessionDate = cellstr(DateStrings);
    Selection   = 1:numel(SessionDate);
else
    if ischar(SessionDate)
        SessionDate = {SessionDate};
    end
    Selection = nan(1,numel(SessionDate));
    for d = 1:numel(SessionDate)
        Selection(d) = strmatch(SessionDate{d}, DateStrings);
        if isempty(Selection(d))
            error('Specified session date ''%s'' was not found in %s!', SessionDate{d}, HistoryFile);
        end
    end
end


%========================== Load spike quality data
SessionParams.ContactData = [];
if strcmpi(HistoryFormat, '.xls')
    [status, sheets] = xlsfinfo(HistoryFile);
    if numel(sheets) >= 2
        [num2,txt2,raw2] =  xlsread(HistoryFile,sheets{2},'');                         % Read data from Excel file sheet 2
        for d = 1:numel(SessionDate)
            DateIndx = strfind(num2(1,:), datenum(SessionDate{d})-datenum('30-Dec-1899'));
            if ~isempty(DateIndx)
                for e = 1:numel(DateIndx)
                    SessionParams(d).ContactData{e} = num2(3:end,DateIndx(e));
                    SessionParams(d).ContactData{e}(isnan(SessionParams(d).ContactData{e})) = 0;
                end
            else
                SessionParams(d).ContactData{1} = 0;
            end
        end
    end
end

%========================== Return data
HeaderStrings = {'Electrode','ML','AP','Depth','GuideLength'};
for h = 1:numel(HeaderStrings)
    ColumnIndx{h} = find(~cellfun(@isempty, strfind(Header, HeaderStrings{h})));
end
for d = 1:numel(Selection)                                                                                      % For each session selected...
    SessionParams(d).Date                   = datestr(C{Selection(d),1});                                   	% Record session date string
    SessionParams(d).DateString             = DateStrings(Selection(d),:);
    SessionParams(d).DateNum                = datenum(DateStrings(Selection(d),:));
    SessionParams(d).DateIndex              = Selection(d);
    SessionParams(d).NoElectrodes           = numel(find(~cellfun(@isnan, C(Selection(d),3:5:end))));        	% How many electrodes were used?
    for e = 1:SessionParams(d).NoElectrodes                                                                     % For each electrode...
        SessionParams(d).ElectrodeID{e}     = C{Selection(d),ColumnIndx{1}(e)};                             	% Get the elctrode identifier
        SessionParams(d).Target{e}          = [C{Selection(d),ColumnIndx{2}(e)},C{Selection(d),ColumnIndx{3}(e)}];    % Get medial-lateral and anetrior-posterior grid hole coordinates
        SessionParams(d).Depth{e}           = C{Selection(d),ColumnIndx{4}(e)};                               	% Get final tip depth
        SessionParams(d).GuideLength{e}     = C{Selection(d),ColumnIndx{5}(e)};                              	% Get the length of guide tube used
    end
end
