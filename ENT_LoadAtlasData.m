function Atlas = ENT_LoadAtlasData

[EN_Root, ~]     = fileparts(mfilename('fullpath'));
AtlasLabelsFile	 = fullfile(EN_Root, 'Atlases','EN_RhesusAtlasLabels.mat');
if ~exist(AtlasLabelsFile, 'file')
    error('ElectroNav file containing atlas labels (''%s'') was not found!', AtlasLabelsFile);
end
load(AtlasLabelsFile);