function Nii = LoadStructureVolumes(StructureDir)

% LoadStructureVolumes

AllStructures = wildcardsearch(StructureDir, '*.nii');
for S = 1:numel(AllStructures)
    [a,b,c] = fileparts(AllStructures{S});
    StructureNames{S} = b;
end
[Selection,ok] = listdlg('ListString',StructureNames,'SelectionMode','multi','PromptString','Select structures:','ListSize',[300, 200]);
if ok==0
    return;
end
H = waitbar(0,sprintf('Loading selected structure %d of %d...',S,numel(Selection)));
Hh = get(findobj(H,'type','axes'),'title');
for S = 1:numel(Selection)
    waitbar((S-1)/numel(Selection),H);
    set(Hh, 'string', sprintf('Loading selected structure %d of %d...',S,numel(Selection)));
    StructNii(S) = load_nii(AllStructures{Selection(S)});
    VoxelDim(S,:) = StructNii(S).hdr.dime.pixdim(2:4);                     % Get voxel dimensions
    VolumeDim(S,:) = size(StructNii(S).img);                               % Get volume dimensions
end
close(H);

%=============== Check volumes



%=============== Add volumes
OriginalLayers = numel(Layer.MRI);
for S = 1:numel(Selection)
    n = S+OriginalLayers;
    Layer.MRI(n).img = double(StructNii(S)img);                                   	% Save image volume
    Layer.MRI(n).VoxelDim = StructNii(S)hdr.dime.pixdim(2:4);                     	% Get voxel size (mm)
    Layer.MRI(n).DimVox = size(StructNii(S)img);                                 	% Get full volume dimensions (voxels)
    Layer.MRI(n).DimMM = Layer.MRI(n).DimVox.*Layer.MRI(n).VoxelDim;                % Convert volume dim to mm
    Layer.MRI(n).OriginVox = StructNii(S)hdr.hist.originator(1:3);               	% Get origin coordinates (voxels)
    Layer.MRI(n).OriginMM = Layer.MRI(n).OriginVox.*Layer.MRI(n).VoxelDim;        	% Convert origin to mm
end
Layer.StructNames = AllStructures{Selection};
set(Layer.InputHandle(2), 'string', Layer.StructNames);
