function [Status] = ENT_WriteSessionParams(HistoryFile, Params)

%========================== ENT_WriteSessionParams.m ======================
% This function writes the provided session parameters to the specified 
% spreadsheet file, for the date specified.
%
% INPUTS:
%       HistoryFile:    full path of spreadsheet (.xls/ .csv) containing recording history
%       Params:         optional string specifying which date to query
%                       parameters for, in DD-MMM-YYYY format. If not provided,
%                       user must select from a list.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy
% ? Copyleft 2014, GNU General Public License
%==========================================================================

%========================== Check inputs
if nargin == 0
    [file, path] = uigetfile({'*.xls;*.csv'},'Select recording history');
    HistoryFile = fullfile(path, file);
end
if ~exist(HistoryFile,'file')
    error('Session history file %s does not exist!', HistoryFile);
end

%========================== Load recording history data
[a,b,HistoryFormat] = fileparts(HistoryFile);
if exist('datetime.m','file')==2                   	%============ MATLAB R2014a and later
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
else                                                %============ MATLAB R2013b and earlier    
    [num,txt,raw] =  xlsread(HistoryFile,1,'');                 % Read data from Excel file
    Headers = txt{1,:};                                       	% Skip row containing column titles
    num(1,:) = [];                                            	% Remove header NaNs
    Dates = num(:,1)+datenum('30-Dec-1899');                 	% Convert Excel dates to Matlab dates
    DateStrings = datestr(Dates);                               
    C = raw(2:end,:);                                           
end

%========================== Get a date input


%========================== Load spike quality data
if strcmpi(HistoryFormat, '.xls')
    [status, sheets] = xlsfinfo(HistoryFile);
    if numel(sheets) >= 2
        [num,txt,raw] =  xlsread(HistoryFile,sheets{2},'');                                     % Read data from Excel file sheet 2
        for d = 1:numel(SessionDate)                                                            
            DateIndx = strfind(num(1,:), datenum(SessionDate{d})-datenum('30-Dec-1899'));       
            for e = 1:numel(DateIndx)
                SessionParams(d).ContactData{e} = num(3:end,DateIndx(e));                       
                SessionParams(d).ContactData{e}(isnan(SessionParams(d).ContactData{e})) = 0;    
            end
        end
    end
end

Electrode(e).ContactData  	= Params(1).ContactData{e};


%========================== Write data
for d = 1:numel(Selection)                                                                                  % For each session selected...
    SessionParams(d).Date               = datestr(C{Selection(d),1});                                       % Record session date string
    SessionParams(d).DateString         = DateStrings(Selection(d),:);
    SessionParams(d).DateIndex       	= Selection(d);
    SessionParams(d).NoElectrodes       = numel(find(~cellfun(@isnan, C(Selection(d),3:5:end))));         	% How many electrodes were used?
    for e = 1:SessionParams(d).NoElectrodes                                                                 % For each electrode...
        SessionParams(d).Target{e}          = [C{Selection(d),3+((e-1)*5)},C{Selection(d),4+((e-1)*5)}];    % Get medial-lateral and anetrior-posterior grid hole coordinates
        SessionParams(d).Depth{e}           = C{Selection(d),5+((e-1)*5)};                                  % Get final tip depth
        SessionParams(d).ElectrodeID{e}     = C{Selection(d),2+((e-1)*5)};                                  % Get the elctrode identifier
        SessionParams(d).GuideLength{e}     = C{Selection(d),6+((e-1)*5)};                                  % Get the length of guide tube used
    end
end
