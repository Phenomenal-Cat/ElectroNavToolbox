% [SelectedIndx, SelectedStructs] = GetStructureIndex(Searchterms, Bilateral)

%============================ GetStructureIndex ===========================
% Reads the atlas structure labels for the NeuroMaps Atlas from text file, 
% and presents a list of structure names for the user to select from. The
% structure index values for selected structures are returned. 
%
% INPUTS (optional):
%   Bilateral:  0 = list each hemispheres structures individually. 
%               1 = structures from both hemispheres are automatically selected.
%   Searchterms: cell containing strings of search terms to include.
%
% EAMPLE:
%   [SelectedIndx, SelectedStructs] = GetStructureIndex('Lateral genic', 1);
%
% REVISION HISTORY:
%   01/10/2013 - Written by APM (murphyap@mail.nih.gov)
%   21/08/2014 - Updated to allow search terms
%   11/10/2014 - Updated to allow flexibility in search term format
%==========================================================================

function [SelectedStructIndx, SelectedStructs] = GetStructureIndex(Searchterms, Bilateral)

if nargin == 0                                                  % If no input variable was provided...
    Bilateral = 1;                                              % Default to bilateral structure selection
    Searchterms = [];                                            
end

%============= Load list of structure names
LabelsFile = 'Atlases/inia19/inia19-NeuroMaps.txt';          	% Specify filename of NeuroMaps structure data
if exist(LabelsFile)==0
    error('Atlas structure file ''%s'' was not found!\n', LabelsFile);
end
Stringformat = '%d %s %d %d %d %d';                             % Specify format of NeuroMaps structure data
fid = fopen(LabelsFile);                                        % Open file
txt = textscan(fid, Stringformat);                              % Read strcuture data from file
fclose(fid);                                                    
StructureNames = sort(txt{2});                                  % Sort structure names alphabetically
if Bilateral == 1                                               % For bilateral structure selection...
    StructureNames = StructureNames(strncmp('l_',txt{2},2));    
    StructureNames = cellfun(@(s) s(3:end), StructureNames, 'UniformOutput', false);
end

%============= Find selected structures
if ~isempty(Searchterms)
    if ischar(Searchterms)
        Searchterm = Searchterms;
        clear Searchterms;
        Searchterms{1} = Searchterm;
    end
    if ~iscell(Searchterms)
        error('Search term input is type: ''%s''. Accepted input types are: cell array or string!\n', class(Searchterms));
    end
    Searchterms = lower(Searchterms);                           % Make search terms lower-case
    Selection = [];
    for t = 1:numel(Searchterms)
      	Searchterms{t}(regexp(Searchterms{t},' ')) = '_';       % Replace any white space in search term with underscore
        IndicesA = strfind(StructureNames,Searchterms{t});      % Find cells containing search term t
        IndicesB = find(not(cellfun('isempty', IndicesA)));     % Find non-empty cells
        Selection = [Selection, IndicesB'];                    	% Add indices to selection list
    end
    
else
    [Selection,ok] = listdlg(   'ListString',StructureNames,...
                                'SelectionMode','multi',...
                                'ListSize',[300 300],...
                                'PromptString','Select structure(s):');
    if ok == 0
        SelectedStructIndx = [];
        return;
    end
end

%============= Return atlas indices for selected structures
si = 1;                                                         % Set structure index count to 1
for s = 1:numel(Selection)                                    	% For each stucture selected...
    SelectedStructs{s} = StructureNames{Selection(s)};        	% Find selected structure name
    Selected = strfind(txt{2},SelectedStructs{s});              % Get index in structure name list
    Index = find(not(cellfun('isempty', Selected)));            
    for i = 1:numel(Index)                                      % In case number of indices is more than 1...
        SelectedStructIndx{s}(i) = txt{1}(Index(i));            
        SelectedIndx(si) = txt{1}(Index(i));                  	% Find NeuroMaps index for selected structure(s)
        si = si+1;                                              % Advance structure count
    end
end
end