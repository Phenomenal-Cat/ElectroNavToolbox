
% RecombineWarpedStructures

StructDir = '/Users/aidanmurphy/Pulvinar_project/ElectroNavToolbox/Subjects/Layla/StructureVolumes';
OutputDir = '/Users/aidanmurphy/Pulvinar_project/ElectroNavToolbox/Subjects/Layla/';
OutputFile = 'warped_Layla_INIA19.nii';
WarpedVols = wildcardsearch(StructDir, 'warped_*.nii');
StructNames = {};
Thresh = 0.5;                   % Threshold for determining whether a voxel belongs to a structure (post-interpolation)
WarpedVols = uipickfiles('FilterSpec',fullfile(StructDir,'*nii'),'Prompt','Select session to plot');

for w = 1:numel(WarpedVols)
    [a,Filename,c] = fileparts(WarpedVols{w});
    StructNames{w} = Filename(8:end);
    TempNii = load_nii(WarpedVols{w});
    if w==1
        Combined = TempNii;                              % Copy nifti structure
    else
        Combined.img(TempNii.img>=Thresh) = w;           % Add new structure to volume
    end
end
save_nii(Combined, fullfile(OutputDir,OutputFile));                     % Save new .nii file
save(fullfile(OutputDir,'warped_Layla_INIA19.mat'),'StructNames');      % Save structure names to .mat file