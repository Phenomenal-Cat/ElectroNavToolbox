
%========================== ENT_AlignCenters.m ============================ 
% This subfunction is used to align the origins of two volumes that are
% already coregistered but have different resolutions, and/or bounding box 
% sizes. For example, if a T1-weighted (MDEFT) and T2-weighted (FLASH) 
% structural volumes were acquired in the same session at different resolutions
% then they should already be spatially registered.
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function ENT_AlignCenters(NiiFile1, NiiFile2, verbose)

if nargin<2
    [File1, Path] = uigetfile({'*.nii;*.img'},'Select first MR volume');
    NiiFile1 = fullfile(Path, File1);
    [File2, Path] = uigetfile({'*.nii;*.img'},'Select second MR volume', Path);
    NiiFile2 = fullfile(Path, File2);
    NiiFiles = {NiiFile1, NiiFile2};
    verbose = 1;
else
    verbose = 0;
end

%===== Load volumes from files
if ~isstruct(NiiFile1)
    Nii1 = load_nii(NiiFile1);
    Nii2 = load_nii(NiiFile2);
end

%===== Check volumes have same resolution
if Nii1.hdr.dime.pixdim(2:4)~= Nii2.hdr.dime.pixdim(2:4)
    NewRes  = min([Nii1.hdr.dime.pixdim(2:4), Nii2.hdr.dime.pixdim(2:4)]);
    Indx    = find([Nii1.hdr.dime.pixdim(2:4), Nii2.hdr.dime.pixdim(2:4)]==NewRes);
    
    NewNiiFile = NiiFiles{Indx};
    reslice_nii(NiiFiles{Indx}, NewNiiFile, NewRes);
    NiiFiles{Indx} = load_nii(NewNiiFile);
    
end

%===== Plot data to figure
if verbose == 1
    FH = figure;
    ax(1) = subplot(1,3,1);
    imagesc(Nii1.img(:,:,Nii1.hdr.hist.originator(3)));
    axis equal tight;
    hold on;
    plot(xlim, repmat(Nii1.hdr.hist.originator(1),[1,2]),'-w');
    plot(repmat(Nii1.hdr.hist.originator(2),[1,2]), ylim,'-w');
    title(NiiFile1);

    ax(2) = subplot(1,3,2);
    imagesc(Nii2.img(:,:,Nii2.hdr.hist.originator(3)));
    axis equal tight;
    hold on;
    plot(xlim, repmat(Nii2.hdr.hist.originator(1),[1,2]),'-w');
    plot(repmat(Nii2.hdr.hist.originator(2),[1,2]), ylim,'-w');
    title(NiiFile2);

    ax(3) = subplot(1,3,3);
    imagesc(Nii1.img(:,:,Nii1.hdr.hist.originator(3)));
    axis equal tight;
    hold on;
end


%===== Calculate volume dimensions
SizeDiff    = Nii1.hdr.dime.dim(2:4)-Nii2.hdr.dime.dim(2:4);
Offset      = Nii1.hdr.hist.originator(1:3)-Nii2.hdr.hist.originator(1:3);

if all(SizeDiff>0)                                      % If volume 1 is larger than volume 2 in all dimensions...
    NewVol = zeros(size(Nii1.img));                     % Create new empty volume
    Xrange = Offset(1)+(1:Nii2.hdr.dime.dim(2));   
    Yrange = Offset(2)+(1:Nii2.hdr.dime.dim(3));
    Zrange = Offset(3)+(1:Nii2.hdr.dime.dim(4));
    NewVol(Xrange, Yrange, Zrange) = Nii2.img;          % Insert original volume
    
    Nii2.img =  NewVol;
    Nii2.hdr.dime.dim(2:4) = size(Nii1.img);
    Nii2.hdr.hist.originator(1:3) = Nii1.hdr.hist.originator(1:3);
    Filename = sprintf('%s_Aligned.nii', NiiFile2(1:strfind(NiiFile2,'.nii')-1));
end
imagesc(Nii2.img(:,:,Nii1.hdr.hist.originator(3)));
title(Filename);

KbWait;
save_nii(Nii2, Filename);                               % Save new volume