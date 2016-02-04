function Defaults = ENT_LoadDefaults(SubjectID, DefaultFilename)

%============================ ENT_LoadDefaults.m ==========================
% This function loads default parameters for a given subject, for the 
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2014-2016, GNU General Public License
%========================================================================== 

if ~exist('DefaultFilename' , 'var')                                                            % If default filename input wasn't provided...
    EN_Root         = fileparts(mfilename('fullpath'));                                         
    [t, CompName]   = system('hostname');                                                       % Get computer ID
    CompName(strfind(CompName, char(10))) = [];                                                 % Remove white space    
    CompName(strfind(CompName, char(32))) = [];
    DefaultFilename = fullfile(EN_Root,'Params', sprintf('ParamsFile_%s.mat', CompName));
end

if exist(DefaultFilename,'file')~=2                                                             % If default filename doesn't exist...
    msg = sprintf('Parameters file named ''%s'' was not found!', DefaultFilename);
    h = msgbox(msg, 'Warning!');
    uiwait(h);
    [file, path] = uigetfile('*.mat','Select parameters file', fullfile(EN_Root,'Params'));     
    DefaultFilename = fullfile(path, file);
end

load(DefaultFilename);                                                                          % Load default parameters
if exist('SubjectID','var') && ~isempty(SubjectID)                                            	% If SubjectID input was provided...
    SubjectIndx = find(~cellfun(@isempty, strfind({Defaults.SubjectID},SubjectID)));            % Search for that ID
end
if ~exist('SubjectID','var') || isempty(SubjectID)
    SubjectIndx = listdlg('ListString', {Defaults.SubjectID}, 'Name','Select subject ID');      % Ask user to select from available IDs
    if isempty(SubjectIndx)
        Defaults = [];
        return
    end
end
if isempty(SubjectIndx)                                                                         % Otherwise...
    msg = sprintf('No entry was found for subject ID ''%s'' in %s!', SubejctID, DefaultFilename);
    h = msgbox(msg, 'Warning!');
    uiwait(h);
    SubjectIndx = listdlg('ListString', {Defaults.SubjectID}, 'Name','Select from available subject IDs');      % Ask user to select from available IDs
  	if isempty(SubjectIndx)
        Defaults = [];
        return
    end
end
Defaults = Defaults(SubjectIndx);                                                               
    

    