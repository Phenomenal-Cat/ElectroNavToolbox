function Slice = ENT_GetSlice(Subject, Dim, SlicePos)

%========================== ENT_GetSlice.m ================================
% This function returns slice images from an MRI volume for the specified
% subject, slice orientation and position.
%
% INPUTS:   Subject:    string specifying subject ID
%           Dim:        dimension of volume to slice - i.e. slice orientation
%           SlicePos:   position of slice in mm relative to volume origin
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2016, GNU General Public License
%==========================================================================

Paths = ENT_LoadDefaults(Subject);
nii = load_nii(Paths.MRI);

Slice.Origin  = nii.hdr.hist.originator(1:3);
Slice.Scale   = nii.hdr.dime.pixdim(2:4);

for n = 1:3
    Indx{n} = 1:size(nii.img, n);
end
Indx{Dim} = Slice.Origin(Dim) + SlicePos/Slice.Scale(Dim);
Slice.Im = squeeze(nii.img(Indx{1},Indx{2},Indx{3}));