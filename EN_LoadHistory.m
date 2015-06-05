%============================= EN_LoadHistory.m ===========================
% This subfunction loads parameters for all previous recording sessions from 
% the spreadsheet file (.xls or .csv format) provided as an input.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy © Copyleft 2015, GNU General Public License
%==========================================================================

function Hist = EN_LoadHistory(HistoryFile)

if nargin ==0 || exist(HistoryFile, 'file') ~=2                 % If recording history file was not found...
    [file, path] = uigetfile({'*.xls;*.csv'}, 'Select recording history');
    HistoryFile = fullfile(path, file);
end

[a,b,HistoryFormat] = fileparts(HistoryFile);                   % Check file format   
if strcmpi(HistoryFormat, '.xls')                               % If file format was Excel...
    if exist('readtable.m','file')                              % For MATLAB R2014a and later...
        T = readtable(HistoryFile);
        T.Date = datetime(T.Date,'ConvertFrom','excel');
        Hist.DateStrings = char(datetime(T.Date,'format','dd-MMM-yyyy'));   
        
    else                                                        % MATLAB R2013b and earlier...    
        [num,txt,raw] =  xlsread(HistoryFile,1,'');             % Read data from Excel file
        Headers = txt{1,:};                                  	% Skip row containing column titles
        num(1,:) = [];                                       	% Remove nans
        Dates = num(:,1)+datenum('30-Dec-1899');             	% Convert Excel dates to Matlab dates
      	Hist.DateStrings = datestr(Dates);                     
    end
    
elseif strcmpi(HistoryFormat, '.csv')                           % If file format was csv...
    fid = fopen(HistoryFile,'rt');
    Headers = textscan(fid, '%s%s%s%s%s%s%s\n', 1, 'delimiter', ',');
    data = textscan(fid, '%f %d %d %f %f %f %s','headerlines',1,'delimiter',',');
    fclose(fid);
    Dates = data{1};   
    Hist.DateStrings = datestr(Dates); 
end           

end