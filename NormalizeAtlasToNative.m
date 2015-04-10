function WarpFile = NormalizeAtlasToNative(NativeFile, Template, LesionMask)

%===================== NormalizeAtlasToNative.m ===========================
% 
%
%
%
%
%==========================================================================

%============ Default settings
if nargin == 0
    SubjectDir = fullfile(cd,'/Subjects/Layla');
    NativeFile = 'Layla_GridScan_ACPC_BET_Masked.nii';
    LesionMask = 'Layla_ACPC_LesionMask.nii';
    Template = 'inia19-t1-brain.nii';
    AtlasFile = 'inia19-NeuroMaps.nii';
end

TemplateSmooth = 0;
SourceSmooth = 8;



%=========== Check contrast similarity between volumes
NativeNii = load_nii(fullfile(SubjectDir,NativeFile));
TemplateNii = load_nii(fullfile(SubjectDir,Template));

figure;
sh(1) = subplot(2,2,1);
hist(double(NativeNii.img(:)),100);
ylabel('Frequency');
title('Native image volume');

sh(2) = subplot(2,2,3);
hist(double(TemplateNii.img(:)),100);
ylabel('Frequency');
xlabel('Voxel intensities');
linkaxes(sh, 'y');
title('Template image volume');
set(sh([1,2]),'ylim',[0 100000]);

sh(3) = subplot(2,2,2);
imagesc(squeeze(NativeNii.img(NativeNii.hdr.hist.originator(1),:,:)));
colormap gray;
colorbar;

sh(4) = subplot(2,2,4);
imagesc(squeeze(TemplateNii.img(TemplateNii.hdr.hist.originator(1),:,:)));
colormap gray;
set(sh([3,4]),'dataaspect',[1 1 1]);
colorbar;

choice = questdlg('Proceed with normalization?','Check volumes','Yes','No','No');
if strcmpi(choice,'No')
    return;
end

%============ CalculateBoundingBox
MRIVoxelSize = NativeNii.hdr.dime.pixdim(2:4);
MRIOrigin = NativeNii.hdr.hist.originator(1:3)*MRIVoxelSize(1);
MRIVolumeSize = size(NativeNii.img)*MRIVoxelSize(1);
BoundingBox(1,:) = -MRIOrigin;
BoundingBox(2,:) = MRIVolumeSize-MRIOrigin


%=========== Call SPM normalization tool
% spm fMRI;



WarpFile = [Template(1:strfind(Template,'.nii')-1), '_sn.mat'];
load(fullfile(SubjectDir,WarpFile));


params = spm_normalise(VG,VF,matname,VWG,VWF,flags)



