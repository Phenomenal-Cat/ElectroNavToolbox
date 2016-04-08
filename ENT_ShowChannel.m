function [XYZ] = ENT_ShowChannel(SubjectID, Date, Ch)

%============================ ENT_ShowChannel.m ===========================
% Plots the location of the requested cell over the available anatomy for
% the specified subject. 
%
%==========================================================================


[ContactCoords, SessionParams] = EN_GetContactCoordinates(Date, SubjectID);
XYZ         = ContactCoords(1,:,Ch);
AxesLims    = [0, 20; -30, 10; -15, 20];
if strcmpi(SubjectID, 'Layla')
    AxesLims(1,:) = -AxesLims(1,[2,1]);
end
FillStruct  = 1;
MaskAlpha   = 0.5;
fh  = figure;
axh = tight_subplot(1,3,0.02,0.02,0.02);
for i = 1:3
    [Anatomy, Structures, Outlines] = EN_GetAnatomySlices(SubjectID, i, XYZ(i), AxesLims, 3, 0);
    axes(axh(i));
    AnatomyRGB	= repmat(double(Anatomy)/double(max(max(Anatomy))), [1,1,3]);
    imh(i)      = image(AnatomyRGB);
    hold on;
    axis xy tight equal
    StructColors = jet(numel(Structures));
    
    %============= Draw structures
    for N = 1:numel(Structures)
        if FillStruct == 1 
            StructBinary{N} = double(Structures{N})/double(max(max(Structures{N})));        	% Normalize mask values
            imsh(i,N) = imagesc(StructBinary{N}*N,'alphadata', StructBinary{N}*MaskAlpha);       	% Plot structure filled
        end
        size(Outlines)
%             for k = 1:size(Outlines, 3)
%                 StructLineH(i,N) = plot(Outlines{1,N,k}(:,1), Outlines{1,N,k}(:,2),'-w','linewidth',1);     % Plot structure boundary outline
%             end
% %         else
% %             StructLineH(i,N) = plot(Outlines{1,N}(:,1), Outlines{1,N}(:,2),'-w','linewidth',1);     % Plot structure boundary outline
% %         end
%     	set(StructLineH(i,N), 'color', StructColors(N,:));
    end
    
  	%============= Set axes limits
    switch i
        case 1
            set(imh(i),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
            set(imsh(i,:),'xdata', AxesLims(2,:),'ydata', AxesLims(3,:));
            ylabel('Z (mm)');
            xlabel('Y (mm)');
        case 2
            set(imh(i),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
            set(imsh(i,:),'xdata', AxesLims(1,:),'ydata', AxesLims(3,:));
            xlabel('X (mm)');
        case 3
            set(imh(i),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
            set(imsh(i,:),'xdata', AxesLims(1,:),'ydata', AxesLims(2,:));
          	ylabel('Y (mm)');
            xlabel('X (mm)');
    end
    
	%============= Draw cross hairs
 	switch i
        case 1
            plot(AxesLims(2,:), repmat(XYZ(3),[1,2]), '-g');
            plot(repmat(XYZ(2),[1,2]), AxesLims(3,:), '-g');
        case 2
            plot(AxesLims(1,:), repmat(XYZ(3),[1,2]), '-g');
            plot(repmat(XYZ(1),[1,2]), AxesLims(3,:), '-g');
        case 3
            plot(AxesLims(1,:), repmat(XYZ(2),[1,2]), '-g');
            plot(repmat(XYZ(1),[1,2]), AxesLims(2,:), '-g');
    end
        
        
end

%============= 
axes(axh(1));
title(sprintf('XYZ = [%.2f, %.2f, %.2f] mm', XYZ(1), XYZ(2), XYZ(3)), 'fontsize', 18, 'horizontalalignment', 'left');

