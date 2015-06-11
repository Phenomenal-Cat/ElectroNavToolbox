function CombineVolumes(VolumeFiles, OutputFile)

if nargin==0
    VolumeDir = 'Subjects/Layla/StructureVolumes/';
    [filenames, pathname, filt] = uigetfile('*.nii', 'Select volumes to render', VolumeDir, 'multiselect','on');
    for f = 1:numel(filenames)
        VolumeFiles{f} = fullfile(pathname, filenames{f});
    end
    [filename, pathname] = uiputfile('*.nii', 'Save combined volume as:', pathname);
    OutputFile = fullfile(pathname, filename);
end

for s = 1:numel(VolumeFiles)
    nii = load_nii(VolumeFiles{s});
    if s==1
        AllNii = nii;
    else
        AllNii.img = AllNii.img+nii.img;
    end
end
save_nii(AllNii, OutputFile);