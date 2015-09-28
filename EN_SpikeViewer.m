
%=========================== EN_SpikeViewer.m =============================
% This GUI function plots low-level information about spikes for a population
% of neurons. Different variables can be plotted against each other and
% individual cells can be interactively selected via the main plot or
% through the GUI menus, in order to view plots of spike parameters for the
% selected cell.
%
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

function EN_SpikeViewer(DataFile)
global Fig Data Waveform Burst SpikeData

[ENroot, b, c]= fileparts(mfilename('fullpath'));
addpath(genpath(ENroot));
if nargin == 0
    Data.File = '/Volumes/PROJECTS/murphya/Physio/Analysis/SpikeAnalysis/SpikeWaveformData_20140130-20150319.mat';
end

%================= LOAD SPIKE WAVEFORM DATA
FH = EN_About(1);
load(Data.File);
Data.Dates = unique({SpikeData.Date});
Data.Channels = strread(num2str(1:24),'%s');
Data.Cells = strread(num2str(1:3),'%s');
for i = 1:numel(Data.Dates)
    Fig.SpikeIdMat(find(ismember({SpikeData.Date},Data.Dates{i})),1) = i;
end
Fig.SpikeIdMat(:,[2,3]) = [cell2mat({SpikeData.Channel})', cell2mat({SpikeData.Cell})'];
Fig.Current.DateNo = 1;
Fig.Current.ChannelNo = 1;
Fig.Current.CellNo = 1;
close(FH);

%===================== OPEN FIGURE WINDOW
Fig.Background = repmat(0.75,[1,3]);                                            % Set figure background color  
Fig.Handles.Figure = figure('Name',sprintf('ElectroNav%c - Spike Viewer',char(169)),... 	% Open a figure window with specified title
                    'Color',Fig.Background,...                                  % Set the figure window background color
                    'Renderer','OpenGL',...                                     % Use OpenGL renderer
                    'units','normalized',...                                    % Set units to normalized
                    'Position', [0 0 1 1],...                                 	% position figure window to fit fullscreen
                    'NumberTitle','off',...                                     % Remove figure number from title
                    'IntegerHandle','off');                                     % Don't use integer handles
%                 'menu','none','toolbar','none',...                          % Remove toolbar etc

%===================== INITIALIZE GUI CONTROLS
Fig.Current.DateNo = 1;
Fig.Current.ChannelNo = 1;
Fig.Current.CellNo = 1;
Fig.Pannel.Handle = uipanel('BackgroundColor',Fig.Background,'Units','normalized','Position',[0.7,0.1,0.25,0.75],'Title','Current cell','FontSize',14);
Fig.Pannel.Labels = {'Data file','Session date','Channel #','Cell #','X-axis','Y-axis'};
Fig.Pannel.InputType = {'pushbutton','popup','popup','popup','popup','popup'};
Fig.AxisLabels = {'Frequency','Spike width (ms)','Spike amplitude (uV)','Median inter-spike interval (ms)','Median inter-burst interval (ms)'};
Fig.Pannel.Inputs = {Data.File, Data.Dates, Data.Channels, Data.Cells, Fig.AxisLabels, Fig.AxisLabels};
set(Fig.Pannel.Handle,'units','pixels')
PannelSize = get(Fig.Pannel.Handle,'position');
Ypos = PannelSize(4)-(0:25:150)-50;
for n = 1:numel(Fig.Pannel.Labels);
    Fig.Pannel.LabelsH(n) = uicontrol('Style','Text','String',Fig.Pannel.Labels{n},'HorizontalAlignment','Left','pos',[10, Ypos(n), 80, 25],'parent',Fig.Pannel.Handle);
    Fig.Pannel.InputsH(n) = uicontrol('Style',Fig.Pannel.InputType{n},'String',Fig.Pannel.Inputs{n},'pos',[100, Ypos(n), 150, 25],'parent',Fig.Pannel.Handle,'callback',{@MenuInput,n});
end
ThreshBoxPos = [270, Ypos(5), 60, 20; 350, Ypos(5), 60, 20; 270, Ypos(6), 60, 20; 350, Ypos(6), 60, 20]; 
for th = 1:size(ThreshBoxPos,1)
    Fig.Pannel.ThreshInput(th) = uicontrol('Style', 'edit', 'String', '0','parent',Fig.Pannel.Handle,'position',ThreshBoxPos(th,:), 'callback',{@MenuInput,n+th});
end
Fig.Xaxis = 2;
Fig.Yaxis = 3;
set(Fig.Pannel.LabelsH, 'backgroundcolor', Fig.Background);
set(Fig.Pannel.InputsH(5), 'value', Fig.Xaxis);
set(Fig.Pannel.InputsH(6), 'value', Fig.Yaxis);

%===================== DRAW INDIVIDUAL CELL AXES
Fig.Pannel.AxPos(1,:) = [0.15, 0.6, 0.7, 0.18];
Fig.Pannel.AxPos(2,:) = [0.15, 0.35, 0.7, 0.18];
Fig.Pannel.AxPos(3,:) = [0.15, 0.1, 0.7, 0.18];
for a = 1:size(Fig.Pannel.AxPos,1)
    Fig.Pannel.AxH(a) = axes('units','normalized','position',Fig.Pannel.AxPos(a,:),'parent',Fig.Pannel.Handle);
end

%===================== DRAW MAIN AXES
Fig.AxPos(1,:) = [0.1, 0.15, 0.55, 0.75];
Fig.AxH = axes('units','normalized','position',Fig.AxPos(1,:));
for i = 1:numel(SpikeData)
    Fig.Data.Handle(i) = plot(Waveform(i).PeakToTroughDuration, abs(Waveform(i).PeakToTroughAmp) ,'.b','ButtonDownFcn',{@DataSelect, i});
    hold on;
end
set(Fig.Data.Handle, 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);
box off;
grid on;
xlabel(Fig.AxisLabels{Fig.Xaxis},'fontsize',18);
ylabel(Fig.AxisLabels{Fig.Yaxis},'fontsize',18);
Fig.Xlims = get(Fig.AxH, 'xlim');
Fig.Ylims = get(Fig.AxH, 'ylim');
for i = 1:2
    set(Fig.Pannel.ThreshInput(i), 'String', num2str(Fig.Xlims(i)));
    set(Fig.Pannel.ThreshInput(i+2), 'String', num2str(Fig.Ylims(i)));
end

PlotCellData;   %===== Plot data for current selected cell

end

%========================================================================= 
function DataSelect(objectHandle, eventData, Indx)
global Fig Waveform Data SpikeData
    Fig.Current.CellIndx = Indx;                                                                                        % Get cell index number
    set(Fig.Data.Handle, 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);                          % Reset all cells
    set(Fig.Data.Handle(Fig.Current.CellIndx),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',30);  	% Highlight selected
    
    %=========== UPDATE UI CONTROLS
    Fig.Current.DateNo = Fig.SpikeIdMat(Fig.Current.CellIndx, 1);               
    Fig.Current.ChannelNo = Fig.SpikeIdMat(Fig.Current.CellIndx, 2);
    Fig.Current.CellNo = Fig.SpikeIdMat(Fig.Current.CellIndx, 3);
  	Data.Channels = num2str(unique(Fig.SpikeIdMat(find(Fig.SpikeIdMat(:,1)==Fig.Current.DateNo),2)));
    Data.Cells = num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3)));
 	set(Fig.Pannel.InputsH(2), 'value', Fig.Current.DateNo);
    set(Fig.Pannel.InputsH(3), 'string', Data.Channels,'value',Fig.Current.ChannelNo);
    set(Fig.Pannel.InputsH(4), 'string', Data.Cells,'value',Fig.Current.CellNo);

    PlotCellData;                   
end

%======================== REPLOT POPULATION DATA ========================== 
function PlotPopData
global Fig Data Waveform SpikeData Burst
    axes(Fig.AxH(1));                                           % Select main plot
    xlabel(Fig.AxisLabels{Fig.Xaxis},'fontsize',18);            % Update axis labels
    ylabel(Fig.AxisLabels{Fig.Yaxis},'fontsize',18);
    
    for i = 1:numel(Fig.Data.Handle)
        switch Fig.Xaxis
            case 2	%================== SPIKE WIDTH
                set(Fig.Data.Handle(i), 'XData', Waveform(i).PeakToTroughDuration);
            case 3	%================== SPIKE AMPLITUDE
                set(Fig.Data.Handle(i), 'XData', abs(Waveform(i).PeakToTroughAmp));
            case 4	%================== ISI
                set(Fig.Data.Handle(i), 'XData', Burst(i).ISIPercentiles(3));
            case 5	%================== IBI
                set(Fig.Data.Handle(i), 'XData', Burst(i).IBIPercentiles(3));  
        end

        switch Fig.Yaxis
            case 2	%================== SPIKE WIDTH
                set(Fig.Data.Handle(i), 'YData', Waveform(i).PeakToTroughDuration);   
            case 3
                set(Fig.Data.Handle(i), 'YData', abs(Waveform(i).PeakToTroughAmp));
            case 4
                set(Fig.Data.Handle(i), 'YData', Burst(i).ISIPercentiles(3));
            case 5
                set(Fig.Data.Handle(i), 'YData', Burst(i).IBIPercentiles(3)); 
        end
    end
    
    %================== HISTOGRAM?
    if Fig.Yaxis == 1
        set(Fig.AxH(1), 'visible','off');
        set(Fig.Data.Handle, 'visible','off');
        [n,x] = hist(cell2mat(get(Fig.Data.Handle, 'XData')), 100);
        if isfield(Fig,'HistAx') && ishandle(Fig.HistAx)
            set(Fig.HistAx, 'visible','on');
         	delete(Fig.HistH);
        else
            Fig.HistAx = axes('units','normalized','position',Fig.AxPos(1,:));
        end
       	Fig.HistH = bar(x,n);
        set(Fig.HistAx, 'box','off','tickdir','out','fontsize',14,'ylim',[0, max(n)]);
      	xlabel(Fig.AxisLabels{Fig.Xaxis},'fontsize',18);                    % Update axis labels
        ylabel(Fig.AxisLabels{Fig.Yaxis},'fontsize',18);
        Fig.Ylims = [0, max(n)];
        Fig.Xlims = [min(x), max(x)];

    else
        set(Fig.AxH(1), 'visible','on');
        set(Fig.Data.Handle, 'visible','on');
        if isfield(Fig,'HistAx') && ishandle(Fig.HistAx)
            set(Fig.HistAx, 'visible','off');
            delete(Fig.HistH);
        end
        Fig.Xlims = [min(cell2mat(get(Fig.Data.Handle, 'xdata'))), max(cell2mat(get(Fig.Data.Handle, 'xdata')))];
        Fig.Ylims = [min(cell2mat(get(Fig.Data.Handle, 'ydata'))), max(cell2mat(get(Fig.Data.Handle, 'ydata')))];
    end
    
   	%=========== Update axes limits
   	set(Fig.AxH, 'xlim',Fig.Xlims, 'ylim', Fig.Ylims);
    for i = 1:2
        set(Fig.Pannel.ThreshInput(i), 'String', sprintf('%.2f', Fig.Xlims(i)));
        set(Fig.Pannel.ThreshInput(i+2), 'String', sprintf('%.2f', Fig.Ylims(i)));
    end
    
end

%===================== PLOT SELECTED CELL DATA
function PlotCellData
global Fig Data Waveform Burst
    if isfield(Fig.Current,'Handles') && ishandle(Fig.Current.Handles(1))
        delete(Fig.Current.Handles);
    end
    Fig.Current.CellIndx = find(ismember(Fig.SpikeIdMat, [Fig.Current.DateNo, Fig.Current.ChannelNo, Fig.Current.CellNo],'rows'));
   	set(Fig.Data.Handle, 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);                          % Reset all cells
    set(Fig.Data.Handle(Fig.Current.CellIndx),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',30);  	% Highlight selected
    uistack(Fig.Data.Handle(Fig.Current.CellIndx), 'top');                                                              % Bring selected to front
    
    axes(Fig.Pannel.AxH(1)); 	%================ plot spike waveform
    [Ha Hb Hc] = shadedplot(Waveform(Fig.Current.CellIndx).Times, Waveform(Fig.Current.CellIndx).Mean-Waveform(Fig.Current.CellIndx).Std, Waveform(Fig.Current.CellIndx).Mean+Waveform(Fig.Current.CellIndx).Std);
    hold on;
    Hd = plot(Waveform(Fig.Current.CellIndx).Times, Waveform(Fig.Current.CellIndx).Mean, '-b', 'linewidth', 2);
    set(Ha(2), 'facecolor', [0.5 0.5 1]);
    set([Hb, Hc], 'visible', 'off');
    ph(1) = plot([0 Waveform(Fig.Current.CellIndx).TroughTime], [Waveform(Fig.Current.CellIndx).TroughAmp,Waveform(Fig.Current.CellIndx).TroughAmp], '--r');
    ph(2) = plot([0 Waveform(Fig.Current.CellIndx).PeakTime], [Waveform(Fig.Current.CellIndx).PeakAmp,Waveform(Fig.Current.CellIndx).PeakAmp], '--r');
    ph(3) = plot([Waveform(Fig.Current.CellIndx).PeakTime,Waveform(Fig.Current.CellIndx).TroughTime],[Waveform(Fig.Current.CellIndx).PeakAmp,Waveform(Fig.Current.CellIndx).PeakAmp],'-r','linewidth',5);
    ph(4) = plot(repmat(Waveform(Fig.Current.CellIndx).TroughTime,[1,2]),[Waveform(Fig.Current.CellIndx).PeakAmp,Waveform(Fig.Current.CellIndx).TroughAmp],'--r');
    set(ph(3), 'ButtonDownFcn', {@ChangeWaveTime,1});
    set(ph(4), 'ButtonDownFcn', {@ChangeWaveTime,2});
    axis tight
    
    axes(Fig.Pannel.AxH(2)); 	%================ plot inter-spike interval distribution
    [N,X] = hist(Burst(Fig.Current.CellIndx).ISIs, 0:20:1000);
    ph(5) = bar(X,N);
    hold on;
    ph(6) = plot(repmat(median(Burst(Fig.Current.CellIndx).ISIs),[1,2]),ylim,'-r');
    ph(7) = plot(repmat(mean(Burst(Fig.Current.CellIndx).ISIs),[1,2]),ylim,'-m');
	axis tight
    legend({'','median','mean'});
    
    axes(Fig.Pannel.AxH(3)); 	%================ plot inter-burst interval distribution
    [N,X] = hist(Burst(Fig.Current.CellIndx).IBIs, 0:20:1000);
	ph(8) = bar(X,N);
    hold on;
    ph(9) = plot(repmat(median(Burst(Fig.Current.CellIndx).IBIs),[1,2]),ylim,'-r');
    ph(10) = plot(repmat(mean(Burst(Fig.Current.CellIndx).IBIs),[1,2]),ylim,'-m');
	axis tight
    
    Xlabels = {'Time (us)','ISI (ms)', 'IBI (ms)'};
    Ylabels = {'Voltage (uV)','No. spikes','No. bursts'};
    for a = 1:size(Fig.Pannel.AxPos,1)
        axes(Fig.Pannel.AxH(a))
        xlabel(Xlabels{a}, 'fontsize', 16);
        ylabel(Ylabels{a}, 'fontsize', 16);
    end
    Fig.Current.Handles = [Ha, Hb, Hc, Hd, ph];
	set(Fig.Pannel.AxH, 'box','off');
end


function ChangeWaveTime(hObj, Event, Indx)
global Fig Data Waveform 
    axes(Fig.Pannel.AxH(1));
    [NewX,NewY] = ginput(1);
    Complete = 0;
    set(gcf,'Pointer','crosshair');
    switch Indx
        case 1
%             while Complete <= 100
%                 NewX = get(gca, 'CurrentPoint');
%                 NewX = NewX(1,1);
                NewIndx = find(Waveform(Fig.Current.CellIndx).Times >= NewX(1,1));
                NewY = Waveform(Fig.Current.CellIndx).Mean(NewIndx(1));
                set(Fig.Current.Handles(7), 'xdata', [0, NewX], 'ydata', repmat(NewY,[1,2]));
                set(Fig.Current.Handles(8), 'xdata', [NewX, Waveform(Fig.Current.CellIndx).TroughTime], 'ydata', [NewY NewY]);
                set(Fig.Current.Handles(9), 'ydata', [NewY, Waveform(Fig.Current.CellIndx).TroughAmp]);
%                 Complete =Complete+1;
%             end
         	Waveform(Fig.Current.CellIndx).PeakTime = NewX;
            Waveform(Fig.Current.CellIndx).PeakAmp = NewY;

        case 2
%             while Complete <= 100
%                 NewX = get(gca, 'CurrentPoint');
%                 NewX = NewX(1,1);
                NewIndx = find(Waveform(Fig.Current.CellIndx).Times >= NewX);
                NewY = Waveform(Fig.Current.CellIndx).Mean(NewIndx(1));
                set(Fig.Current.Handles(6), 'xdata', [0, NewX], 'ydata', repmat(NewY,[1,2]));
                set(Fig.Current.Handles(7), 'xdata', [0, NewX]);
                set(Fig.Current.Handles(8), 'xdata', [Waveform(Fig.Current.CellIndx).PeakTime, NewX]);
                set(Fig.Current.Handles(9), 'xdata', [NewX, NewX]);
%                 Complete =Complete+1;
%             end
            Waveform(Fig.Current.CellIndx).TroughTime = NewX;
            Waveform(Fig.Current.CellIndx).TroughAmp = NewY;
    end
    Waveform(Fig.Current.CellIndx).PeakToTroughDuration = Waveform(Fig.Current.CellIndx).TroughTime-Waveform(Fig.Current.CellIndx).PeakTime;
    Waveform(Fig.Current.CellIndx).PeakToTroughAmp = Waveform(Fig.Current.CellIndx).PeakAmp-Waveform(Fig.Current.CellIndx).TroughAmp;           % < IS SIGN CORRECT?
    PlotPopData;
    set(gcf,'Pointer','arrow');
end

%============================= MENU CONTROLS ==============================
function MenuInput(hObj, Event, Indx)
global Fig Data SpikeData
    switch Indx
        case 1  %===================================== SELECT NEW DATA FILE
            [file, path] = uigetfile('*.mat', 'Select spike data file');
            if file ~= 0
                %=========== LOAD NEW DATA
                PreviousDate = Data.Dates(Fig.Current.DateNo);
                load(fullfile(path,file));
                Data.Dates = unique({SpikeData.Date});
                Fig.Current.DateNo = find(ismember(PreviousDate, Data.Dates));
                for i = 1:numel(Data.Dates)
                    Fig.SpikeIdMat(find(ismember({SpikeData.Date},Data.Dates{i})),1) = i;
                end
                Fig.SpikeIdMat(:,[2,3]) = [cell2mat({SpikeData.Channel})', cell2mat({SpikeData.Cell})'];
                if isempty(Fig.Current.DateNo)
                    Fig.Current.DateNo = 1;
                    Fig.Current.ChannelNo = 1;
                    Fig.Current.CellNo = 1;
                end
                
                %=========== UPDATE UI CONTROLS
                Data.Channels = num2str(unique(Fig.SpikeIdMat(find(Fig.SpikeIdMat(:,1)==Fig.Current.DateNo),2)));
                Data.Cells = num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3)));
                set(Fig.Pannel.InputsH(1), 'string', file);
                set(Fig.Pannel.InputsH(2), 'string', Data.Dates,'value',Fig.Current.DateNo);
                set(Fig.Pannel.InputsH(3), 'string', Data.Channels,'value',Fig.Current.ChannelNo);
                set(Fig.Pannel.InputsH(4), 'string', Data.Cells,'value',Fig.Current.CellNo);
              	
                %=========== UPDATE PLOTS
                PlotCellData;
                PlotPopData;
                
            end
            
        case 2  %===================================== SELECT NEW DATE
         	Fig.Current.DateNo = get(hObj,'value');
            Fig.Current.ChannelNo = 1;
            Fig.Current.CellNo = 1;
            
            Data.Channels = strread(num2str(unique(Fig.SpikeIdMat(find(Fig.SpikeIdMat(:,1)==Fig.Current.DateNo),2))'),'%s');
          	Data.Cells = strread(num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3))'),'%s');
            set(Fig.Pannel.InputsH(3), 'string', Data.Channels,'value',Fig.Current.ChannelNo);
          	set(Fig.Pannel.InputsH(4), 'string', Data.Cells,'value',Fig.Current.CellNo);
            
%             Fig.Current.DateIndx = find(ismember({SpikeData.Date}, Data.Dates{Fig.Current.DateNo}));
%             Fig.Current.SpikeIndx = find(ismember(Fig.SpikeIdMat(Fig.Current.DateIndx,:),Fig.SpikeIdMat(:,[2,3]),'rows'));
            
          	set(Fig.Pannel.InputsH(3), 'string', Data.Channels,'value',Fig.Current.ChannelNo);
          	set(Fig.Pannel.InputsH(4), 'string', Data.Cells,'value',Fig.Current.CellNo);
           	PlotCellData;
            
        case 3  %===================================== SELECT NEW CHANNEL
            Fig.Current.ChannelNo = get(hObj,'value');
            Fig.Current.CellNo = 1;
            Data.Cells = strread(num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3))'),'%s');
            set(Fig.Pannel.InputsH(4),'string', Data.Cells,'value',Fig.Current.CellNo);
           	PlotCellData;
            
        case 4  %===================================== SELECT NEW CELL
            Fig.Current.CellNo = get(hObj,'value');
            PlotCellData;
            
        case 5
            Fig.Xaxis = get(hObj,'value');
            PlotPopData;
            
        case 6
            Fig.Yaxis = get(hObj,'value');
            PlotPopData;
            
        case {7,8}
            Fig.Xlims = str2double(get(Fig.Pannel.ThreshInput([1,2]),'string'))';
            set(Fig.AxH, 'xlim', Fig.Xlims);
            
        case {9,10}
            Fig.Ylims = str2double(get(Fig.Pannel.ThreshInput([3,4]),'string'))';
            set(Fig.AxH, 'ylim', Fig.Ylims);
    end
end
