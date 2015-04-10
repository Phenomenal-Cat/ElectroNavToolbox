%========================= GetGridParams.m ================================
% This function returns a structure ('Grid') containing information about
% the physical properties of the recording grid specified by the input
% variable 'GridID'. These parameters are hard coded in this function,
% based on currently available grid formats from Crist Instruments (see
% http://www.cristinstrument.com/products/implant-intro/grids). To add 
% other formats or to customize the defaults for existing formats, edit the 
% parameters listed.
%
% GridID:
%   '19mm_cylindrical':     based on Crist Instrument grid #6-YGD-J1A
%   '35mm_rectangular':     based on custom NIMH/ NEI double grid design
%
% © ElectroNav Toolbox - Aidan Murphy, 2014
%==========================================================================

function Grid = GetGridParams(GridID)

GridFormats = {'19mm_cylindrical','35mm_rectangular'};
[RootDir mfile] = fileparts(mfilename('fullpath'));
LineNumber = 32;                                                            % Line number to edit from
if nargin == 0                                                              % If no input variable...
    Grid = GridFormats;                                                     % Return list of accepted inputs
    return;                                                                 % Exit function
end
Grid.ID = GridID;                                                           % Save grid ID to structure

switch Grid.ID
    
    case '19mm_cylindrical'                 %==================== 19mm cylindrical grid, 17 holes in diameter, 1mm spacing between holes               
        Grid.HoleDiameter = 0.5;                                                    % Grid hole diameter (mm)
        Grid.InterHoleSpacing = 1;                                                  % Distance between centres of adjacent holes (mm)
        Grid.HolesPerColumn = [5 9 11 13 15 15 17 17 17 17 17 15 15 13 11 9 5];     % Number of holes per column
        Grid.HolesPerDim = numel(Grid.HolesPerColumn);                              
        Grid.TotalHoles = sum(Grid.HolesPerColumn);                                 % Total number of grid holes
        Grid.OuterRadius = Grid.InterHoleSpacing*(Grid.HolesPerDim+1)/2;            % Outer radius of grid
        Grid.Width = 10;  
        Grid.Height = 10;
        Grid.GuideTop = 10;
        Grid.RGB = [1 1 0];
        Grid.StlFile = 'Grid1.stl';
        Grid.NiiFile = 'Grid1.nii';
        
    case '35mm_rectangular'              	%==================== 19mm regtangular grid, 17 x 17 holes, 1mm spacing between holes    
        Grid.HoleDiameter = 0.5;                                                    % Grid hole diameter (mm)
        Grid.InterHoleSpacing = 1;                                                  % Distance between centres of adjacent holes (mm)
        Grid.HolesPerColumn = repmat(20,[1,20]);                                    % Number of holes per column
        Grid.HolesPerDim = numel(Grid.HolesPerColumn);                              
        Grid.TotalHoles = sum(Grid.HolesPerColumn);                                 % Total number of grid holes
        Grid.OuterRadius = Grid.InterHoleSpacing*(Grid.HolesPerDim+1)/2;            % Outer radius of grid
        Grid.Width = 10;  
        Grid.Height = 10;
        Grid.GuideTop = 10;
        Grid.RGB = [1 1 0];
        Grid.StlFile = 'Grid2.stl';
        Grid.NiiFile = 'Grid2.nii';

    otherwise
    	Msg = sprintf([ 'The grid ID ''%s'' is not a recognized ID format!\n', ...
                        'Please edit <a href="matlab: opentoline(which('' %s.m''),%d)">%s.m</a> to add parameters for a new grid format.'], Grid.ID, mfilename, LineNumber, mfilename);
        disp(Msg);
        h = warndlg(Msg,'Grid Error');
        return;
        
end

%====================== 3D GRID SETTINGS
Grid.StlFullFile = fullfile(RootDir, 'Grid', Grid.StlFile);        	% Get full path for grid .stl file
[v, f, n, c, stltitle] = stlread(Grid.StlFullFile);                 % Read in .stl file as mesh surface
Grid.faces = f;                                                     
Grid.vertices = v;                                                  
Grid.RGB = [1 1 0];                                                 % Set grid color                                           