
% CompareAtlases

AtlasDir = '/Volumes/APM_1/ElectroNavToolbox/Atlases';
AtlasFiles = {  'NeuroMaps/inia19-NeuroMaps.nii', ...
                'Frey/paxinos_resampled-MNI-full.nii'};
AxisColor = '-w';
AtlasName = {'NeuroMaps','Paxinos'};

figure;
ah = tight_subplot(numel(AtlasFiles), 3, 0.01, 0.01, 0.01);
for A = 1:numel(AtlasFiles)
    Filename = fullfile(AtlasDir, AtlasFiles{A});
    nii(A) = load_nii(Filename);
    if A==1
%         nii(A).img(nii(A).img>1000) = nii(A).img(nii(A).img>1000)-1000;
    elseif A==2
        nii(A).img(nii(A).img>255.4) = 0;
    end
    Origin = nii(A).hdr.hist.originator(1:3);
    VoxelSize = nii(A).hdr.dime.pixdim(2:4);
    
    for V = 1:3
        axes(ah( (A-1)*3 + V));
        [a b c] = size(nii(A).img);
        sX = 1:a;
        sY = 1:b;
        sZ = 1:c;
        if V==1
            sX = Origin(V);
        elseif V==2
            sY = Origin(V);
        elseif V==3
            sZ = Origin(V);
        end
        Slice = rot90(squeeze(nii(A).img(sX,sY,sZ)));
        imagesc(Slice);
        hold on;
        if V==1
            plot(xlim, repmat(Origin(3),[1,2]), AxisColor);
            plot(repmat(Origin(2),[1,2]), ylim, AxisColor);
            xlabel(AtlasName{A},'fontsize',18);
        elseif V==2
          	plot(xlim, repmat(Origin(3),[1,2]), AxisColor);
            plot(repmat(Origin(1),[1,2]), ylim, AxisColor);
        elseif V==3
            plot(xlim, repmat(b-Origin(2),[1,2]), AxisColor);
            plot(repmat(Origin(1),[1,2]), ylim, AxisColor);
        end
        axis equal tight off;
    end
end
linkaxes(ah([1,2]),'y');
linkaxes(ah([4,5]),'y');
Background = [0.5 0.5 0.5];
cm = [Background; jet];
colormap(cm);
set(gcf, 'color', Background);

% export_fig('AtlasComparison.png','-png');
