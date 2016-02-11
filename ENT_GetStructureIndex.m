% AtlasStruct =  ENT_GetStructureIndex(AtlasName, Searchterms, Bilateral)

%======================== ENT_GetStructureIndex.m =========================
% Reads the atlas structure labels for the NeuroMaps Atlas from text file, 
% and presents a list of structure names for the user to select from. The
% structure index values for selected structures are returned. 
%
% INPUTS (optional):
%   AtlasName:      options: 'NeuroMaps','Paxinos','Saleem-Logo'
%   Searchterms:    cell containing strings of search terms to include.
%
% OUTPUT:
%   AtlasStruct.Names:      structure names
%   AtlasStruct.Indices:    numerical value of voxels belonging to structure
%   AtlasStruct.RGB:        Default RGB triplet used for plotting the structure
%
% EAMPLE:
%   AtlasStruct = ENT_GetStructureIndex('neuromaps', 'Lateral genic');
%
% REVISION HISTORY:
%   01/10/2013 - Written by APM (murphyap@mail.nih.gov)
%   21/08/2014 - Updated to allow search terms
%   11/10/2014 - Updated to allow flexibility in search term format
%   07/10/2015 - Updated to return structure names when index is input
%   02/09/2016 - Rewritten to handle multiple atlases
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function [AtlasStruct] = ENT_GetStructureIndex(AtlasName, Searchterms)


%============= Load list of structure names
Atlas = ENT_LoadAtlasData;

if ~exist('AtlasName','var') || ~any(~cellfun(@isempty, strfind(lower({Atlas.name}), lower(AtlasName))))
    Selection = listdlg('ListString',{Atlas.name},'SelectionMode','single','Listsize',[150, 100],'PromptString','Select atlas');
    AtlasName = Atlas(Selection).name; 
end
AtlasIndx = find(~cellfun(@isempty, strfind(lower({Atlas.name}), lower(AtlasName))));   	% Find index of requested atlas

if ~isempty(Atlas(AtlasIndx).StructNames)                                   % If stucture names field contains data
    StructureNames = Atlas(AtlasIndx).StructNames;                          % Use full structure names
elseif isempty(Atlas(AtlasIndx).StructNames)                                % If stucture names field is empty...
    StructureNames = Atlas(AtlasIndx).StructAbbrev;                         % Use structure abbreviations
    if ~isempty(Searchterms)
        fprintf('Warning: Full structure names are not available for the %s atlas. Please specify scructures by their abbreviations.', AtlasName);
    end
end


%============= Find selected structures
if exist('Searchterms','var')
    switch class(Searchterms)
        case 'char'                                                    	% Search term is character...
            Searchterms = {Searchterms};                               	% Put string in cell  
        case 'double'
            Selection   = Searchterms;                          
        otherwise
            if ~iscell(Searchterms)
                error('Search term input is type: ''%s''. Accepted input types are: cell, string or double!\n', class(Searchterms));
            end
    end
    Searchterms	= lower(Searchterms);                                   % Make search terms lower-case
    Selection   = [];                                                   % Initialize variable for selection
    for t = 1:numel(Searchterms)
        Searchterms{t}(regexp(Searchterms{t},' ')) = '_';               % Replace any white space in search term with underscore
        IndicesA    = strfind(StructureNames,Searchterms{t});           % Find cells containing search term t
        IndicesB    = find(not(cellfun('isempty', IndicesA)));          % Find non-empty cells
        Selection   = [Selection, IndicesB'];                         	% Add indices to selection list
    end

elseif ~exist('Searchterms','var')
    StructNamesAlphabet = sortrows(StructureNames);
    [SelectAlpha,ok] = listdlg( 'ListString', StructNamesAlphabet,...
                                'SelectionMode','multi',...
                                'ListSize',[300 300],...
                                'PromptString','Select structure(s):');
  	for s = 1:numel(SelectAlpha)
        Selection(s) = find(strcmp(StructNamesAlphabet{SelectAlpha(s)}, StructureNames));
    end                     
    if ok == 0
        StructIndx = [];
        return;
    end
end

%============= Return names and indices for selected atlas structures
AtlasStruct.Names   = StructureNames(Selection);                        % Find selected structure names
AtlasStruct.RGB     = Atlas(AtlasIndx).StructRGB(Selection, :);         % Get default RGB colors for selected structures
if ~iscell(Atlas(AtlasIndx).StructIndices)
    AtlasStruct.Indices = Atlas(AtlasIndx).StructIndices(Selection);    % Get selected structure indices
elseif iscell(Atlas(AtlasIndx).StructIndices)                           
    for s = 1:numel(Selection)                                          % For each stucture selected...
        AtlasStruct.Indices{s} = Atlas(AtlasIndx).StructIndices{Selection(s)};           
    end
end

end