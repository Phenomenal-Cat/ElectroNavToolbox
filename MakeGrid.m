function MakeGrid()

%============================ MakeGrid.m ==================================
% Generate both a surface mesh (.stl) and volume (.nii) of a recording
% chamber grid based on the specified parameters.
%
% INPUTS:
%   GridShape:     	'round'/'square'
%   HoleDiamater:   hole diameter (mm)
%   HoleSpacing:    inter-hole spacing (mm)
%   HolesPerSide:   [x y] number of holes
%   GrooveType:     type of groove cut into the grid to allow the correct
%                   orientation to be maintained across sessions.
% REQUIREMENTS:
%   - NiftiToolbox 
%   - GPToolbox 
%
%==========================================================================

