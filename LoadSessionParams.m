function [SessionParams] = LoadSessionParams(HistoryFile, SessionDate)

%========================== LoadSessionParams.m ===========================
% This function retreives the session parameters stored in the provided
% spreadsheet file, for the date specified.
%
% INPUTS:
%       HistoryFile:    full path of spreadsheet (.xls/ .csv) containing recording history
%       SessionDate:    optional string specifying which date to query
%                       parameters for, in DD-MMM-YYYY format. If not provided,
%                       user must select from a list.
%
%========================================================================== 

%========================== Check inputs
if nargin == 0
    HistoryFile = 'LaylaElectrodeLocations.xls';
end
if ~exist(HistoryFile,'file')
    error('Session history file %s does not exist!', HistoryFile);
end

%========================== Load recording history data
[a,b,HistoryFormat] = fileparts(HistoryFile);
if exist('readtable.m','file')                   	%============ MATLAB R2014a and later
    if strcmpi(HistoryFormat, '.xls')
        T = readtable(HistoryFile);
        T.Date = datetime(T.Date,'ConvertFrom','excel');
      	DateStrings = datestr(T.Date);
        C = table2cell(T);
        
    elseif strcmpi(HistoryFormat, '.csv') 
%         formatSpec = '%{dd-MMM-yyyy}D%f%f%f%f%f%C';
%         T = readtable(HistoryFile,'Delimiter',',','Format',formatSpec);
     	fid = fopen(HistoryFile,'rt');
        Headers = textscan(fid, '%s%s%s%s%s%s%s\n', 1, 'delimiter', ',');
        data = textscan(fid, '%f %f %f %f %f %f %s','headerlines',1,'delimiter',',');
        fclose(fid);
        Dates = data{1};   
        DateStrings = datestr(Dates);
        C = num2cell(cell2mat(data(1:6)));
        C(:,7) = data{7};
    else
        error('File ''%s'' is not a valid spreadsheet format!', HistoryFile);
    end
else                                                %============ MATLAB R201bb and earlier    
    [num,txt,raw] =  xlsread(HistoryFile,1,'');                 % Read data from Excel file
    Headers = txt{1,:};                                       	% Skip row containing column titles
    num(1,:) = [];                                            	% Remove nans
    Dates = num(:,1)+datenum('30-Dec-1899');                 	% Convert Excel dates to Matlab dates
    DateStrings = datestr(Dates);                        
    C = raw(2:end,:); 
end

%========================== Get a date input
if ~exist('SessionDate','var')
    [Selection,ok] = listdlg('ListString',DateStrings,'SelectionMode','single','PromptString','Select session:');
    if ok==0
        SessionParams = [];
        return;
    end
else
    Selection = strmatch(SessionDate, DateStrings);
    if isempty(Selection)
        error('Specified session date ''%s'' was not found in %s!', SessionDate, HistoryFile);
    end
end

%========================== Return data

SessionParams.DateIndex = Selection;
SessionParams.Target = [C{Selection,2},C{Selection,3}];
SessionParams.Depth = C{Selection,4};
SessionParams.Date = datestr(C{Selection,1});
SessionParams.ElectrodeID = C{Selection,7};
SessionParams.GuideLength = C{Selection, 5};
