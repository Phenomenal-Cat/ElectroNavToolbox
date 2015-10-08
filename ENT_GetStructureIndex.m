% [StructIndx, StructNames, StructRGB] =  ENT_GetStructureIndex(Searchterms, Bilateral)

%======================== ENT_GetStructureIndex.m =========================
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
%   [StructIndx, StructNames, StructRGB] = GetStructureIndex('Lateral genic', 1);
%
% REVISION HISTORY:
%   01/10/2013 - Written by APM (murphyap@mail.nih.gov)
%   21/08/2014 - Updated to allow search terms
%   11/10/2014 - Updated to allow flexibility in search term format
%   07/10/2015 - Updated to return structure names when index is input
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function [StructIndx, StructNames, StructRGB] = ENT_GetStructureIndex(Searchterms, Bilateral)

if nargin == 0                                                  % If no input variable was provided...
    Searchterms = [];                                            
end
if nargin < 2
    Bilateral = 1;                                              % Default to bilateral structure selection
end

%============= Load list of structure names
LabelsFile = 'Atlases/inia19/inia19-NeuroMaps.txt';          	% Specify filename of NeuroMaps structure data
if exist(LabelsFile,'file')==0
    uiwait(msgbox(sprintf('Atlas structure file ''%s'' was not found!\n', LabelsFile),'Error', 'error'));
    [file, path] = uigetfile('.txt', 'Select atlas label file');
    LabelsFile = fullfile(path, file);
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
RGB = [txt{3}, txt{4}, txt{5}];

%============= Find selected structures
ReturnIndx = 1;
if ~isempty(Searchterms)
    switch class(Searchterms)
        case 'char'                                           	% Search term is character...
            Searchterms = {Searchterms};                      	% Put string in cell  
        case 'double'
            ReturnIndx = 0;
            Selection = Searchterms;
        otherwise
            if ~iscell(Searchterms)
                error('Search term input is type: ''%s''. Accepted input types are: cell, string or double!\n', class(Searchterms));
            end
    end
    if ReturnIndx == 1
        Searchterms = lower(Searchterms);                           % Make search terms lower-case
        Selection = [];                                             % Initialize variable for selection
        for t = 1:numel(Searchterms)
            Searchterms{t}(regexp(Searchterms{t},' ')) = '_';       % Replace any white space in search term with underscore
            IndicesA = strfind(StructureNames,Searchterms{t});      % Find cells containing search term t
            IndicesB = find(not(cellfun('isempty', IndicesA)));     % Find non-empty cells
            Selection = [Selection, IndicesB'];                    	% Add indices to selection list
        end
    end
elseif isempty(Searchterms)
    [Selection,ok] = listdlg(   'ListString',StructureNames,...
                                'SelectionMode','multi',...
                                'ListSize',[300 300],...
                                'PromptString','Select structure(s):');
    if ok == 0
        StructIndx = [];
        return;
    end
end


if ReturnIndx == 1  %============= Return atlas indices for selected structures
    si = 1;                                                         % Set structure index count to 1
    for s = 1:numel(Selection)                                    	% For each stucture selected...
        StructNames{s} = StructureNames{Selection(s)};              % Find selected structure name
        Selected = strfind(txt{2},StructNames{s});                  % Get index in structure name list
        Index = find(not(cellfun('isempty', Selected)));            
        for i = 1:numel(Index)                                      % In case number of indices is more than 1...
            StructIndx{s}(i) = txt{1}(Index(i));            
            SelectedIndx(si) = txt{1}(Index(i));                  	% Find NeuroMaps index for selected structure(s)
            StructRGB{s} = RGB(Index(i), :);                        % Find default colormap for selected structures
            si = si+1;                                              % Advance structure count
        end
    end
elseif ReturnIndx == 0	%============= Return structure names for selected indices                                          
    for s = 1:numel(Selection)                                  	% For each index provided...
        StructIndx{s} = Selection(s);                               % Find row containing that index
        StructNames{s} = txt{2}{Selection(s)};                      % Find selected structure name
        StructRGB{s} = RGB(Selection(s), :);                       	% Find default colormap for selected structures
    end
end

end