
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
global Fig Data Waveform Burst

if nargin == 0
    Data.File = 'SpikeWaveformData_20140130-20150319.mat';
else
    
end
load(Data.File);
Data.Dates = unique({SpikeData.Date});
Data.Channels = strread(num2str(1:24),'%s');
Data.Cells = strread(num2str(1:3),'%s');


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
Fig.Pannel.Handle = uipanel('BackgroundColor',Fig.Background,'Units','normalized','Position',[0.7,0.05,0.25,0.75],'Title','Current cell','FontSize',14);
Fig.Pannel.Labels = {'Data file','Session date','Channel #','Cell #'};
Fig.Pannel.InputType = {'pushbutton','popup','popup','popup'};
Fig.Pannel.Inputs = {Data.File, Data.Dates, Data.Channels, Data.Cells};
set(Fig.Pannel.Handle,'units','pixels')
PannelSize = get(Fig.Pannel.Handle,'position');
Ypos = PannelSize(4)-(0:25:100)-50;
for n = 1:numel(Fig.Pannel.Labels);
    Fig.Pannel.LabelsH(n) = uicontrol('Style','Text','String',Fig.Pannel.Labels{n},'HorizontalAlignment','Left','pos',[10, Ypos(n), 80, 25],'parent',Fig.Pannel.Handle);
    Fig.Pannel.InputsH(n) = uicontrol('Style',Fig.Pannel.InputType{n},'String',Fig.Pannel.Inputs{n},'pos',[100, Ypos(n), 150, 25],'parent',Fig.Pannel.Handle,'callback',{@MenuInput,n});
end
set(Fig.Pannel.LabelsH, 'backgroundcolor', Fig.Background);

%===================== DRAW AXES
Fig.Pannel.AxPos(1,:) = [0.15, 0.6, 0.7, 0.15];
Fig.Pannel.AxPos(2,:) = [0.15, 0.35, 0.7, 0.15];
Fig.Pannel.AxPos(3,:) = [0.15, 0.1, 0.7, 0.15];
Xlabels = {'Time (us)','ISI (ms)', 'IBI (ms)'};
Ylabels = {'Voltage (uV)','No. spikes','No. bursts'};
for a = 1:size(Fig.Pannel.AxPos,1)
    Fig.Pannel.AxH(a) = axes('units','normalized','position',Fig.Pannel.AxPos(a,:),'parent',Fig.Pannel.Handle);
  	xlabel(Xlabels{a});
    ylabel(Ylabels{a});
end

Fig.AxPos(1,:) = [0.05, 0.6, 0.25, 0.35];
Fig.AxPos(2,:) = [0.35, 0.6, 0.25, 0.35];
Fig.AxPos(3,:) = [0.05, 0.1, 0.25, 0.35];
Fig.AxPos(4,:) = [0.35, 0.1, 0.25, 0.35];
Xlabels = {'P2P duration (us)','ISI (ms)', '',''};
Ylabels = {'Voltage (uV)','No. spikes','Frequency','Frequency'};
for a = 1:size(Fig.AxPos,1)
    Fig.AxH(a) = axes('units','normalized','position',Fig.AxPos(a,:));
    xlabel(Xlabels{a});
    ylabel(Ylabels{a});
end





PlotPopData;

PlotCellData;


end

%===================== REPLOT POPULATION DATA
function PlotPopData
global Fig Data Waveform
    axes(Fig.AxH(1)); 	%================ plot peak-to-peak duration

    for i = 1:numel(SpikeData)
        Fig.Data.Handle(i) = plot(Waveform(i).PeakToTroughDuration, Waveform(i).ISIs ,'.b','ButtonDownFcn',{@DataSelect, i});
        hold on;
    end
    

end

%===================== PLOT SELECTED CELL DATA
function PlotCellData
global Fig Data Waveform Burst
    if isfield(Fig.Current,'Handles') && ishandle(Fig.Current.Handles)
        delete(Fig.Current.Handles);
    end
    Fig.Current.CellIndx = find(ismember(Fig.SpikeIdMat, [Fig.Current.DateNo, Fig.Current.ChannelNo, Fig.Current.CellNo],'rows'));
    
    
    axes(Fig.Pannel.AxH(1)); 	%================ plot spike waveform
    [Ha Hb Hc] = shadedplot(Waveform.Times, Waveform.Mean-Waveform.Std, Waveform.Mean+Waveform.Std);
    hold on;
    Hd = plot(Waveform.Times, Waveform.Mean, ['-',ColorsLine{CIndx}], 'linewidth', 2);
    set(Ha(2), 'facecolor', ColorsRGB(CIndx,:));
    set([Hb, Hc], 'visible', 'off');
    ph(1) = plot([0 Waveform.TroughTime], [Waveform.TroughAmp,Waveform.TroughAmp], '--r');
    ph(2) = plot([0 Waveform.PeakTime], [Waveform.PeakAmp,Waveform.PeakAmp], '--r');
    ph(3) = plot([Waveform.PeakTime,Waveform.TroughTime],[Waveform.PeakAmp,Waveform.PeakAmp],'-r','linewidth',5);
    ph(4) = plot([Waveform.TroughTime,Waveform.TroughTime],[Waveform.PeakAmp,Waveform.TroughAmp],'--r');
    axis tight
    set(gca, 'box','off');
    
    axes(Fig.Pannel.AxH(2)); 	%================ plot inter-spike interval distribution
    hist(Burst.ISIs, 0:20:1000);
    hold on;
    ph(5) = plot(repmat(median(Burst.ISIs),[1,2]),ylim,'-r');
    ph(6) = plot(repmat(mean(Burst.ISIs),[1,2]),ylim,'-m');
	axis tight
    
    axes(Fig.Pannel.AxH(3)); 	%================ plot inter-burst interval distribution
    hist(Burst.IBIs, 0:20:1000);
    hold on;
    ph(7) = plot(repmat(median(Burst.IBIs),[1,2]),ylim,'-r');
    ph(8) = plot(repmat(mean(Burst.IBIs),[1,2]),ylim,'-m');
	axis tight
    
    Fig.Current.Handles = [Ha, Hb, Hd, ph];
end

%===================== MENU CONTROLS
function MenuInput(hObj, Event, Indx)
global Fig Data
    switch Indx
        case 1
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
                Data.Channels = strread(num2str(unique(Fig.SpikeIdMat(find(Fig.SpikeIdMat(:,1)==Fig.Current.DateNo),2))),'%s');
                Data.Cells = strread(num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3))),'%s');
                set(Fig.Pannel.InputsH(1), 'string', file);
                set(Fig.Pannel.InputsH(2), 'string', Data.Dates,'value',Fig.Current.DateNo);
                set(Fig.Pannel.InputsH(3), 'string', Data.Channels,'value',Fig.Current.ChannelNo);
                set(Fig.Pannel.InputsH(4), 'string', Data.Cells,'value',Fig.Current.CellNo);
              	
                %=========== UPDATE PLOTS
                PlotCellData;
                PlotPopData;
                
            end
        case 2
         	Fig.Current.DateNo = get(hObj,'value');
            Fig.Current.ChannelNo = 1;
            Fig.Current.CellNo = 1;
            
            Data.Channels = strread(num2str(unique(Fig.SpikeIdMat(find(Fig.SpikeIdMat(:,1)==Fig.Current.DateNo),2))'),'%s');
          	Data.Cells = strread(num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3))'),'%s');
            set(Fig.Pannel.InputsH(3), 'string', Data.Channels,'value',Fig.Current.ChannelNo);
          	set(Fig.Pannel.InputsH(4), 'string', Data.Cells,'value',Fig.Current.CellNo);
            
            Fig.Current.DateIndx = find(ismember({SpikeData.Date},Fig.UniqueDates{Fig.Current.DateNo}));
            Fig.Current.SpikeIndx = find(ismember(Fig.SpikeIdMat(Fig.Current.DateIndx,:),Fig.SpikeIdMat(:,[2,3]),'rows'));
            
            set(Fig.Pannel.InputsH([3,4]),'value',1);
           	PlotCellData;
            
        case 3
            Fig.Current.ChannelNo = get(hObj,'value');
            Fig.Current.CellNo = 1;
            Data.Cells = strread(num2str(unique(Fig.SpikeIdMat(find(ismember(Fig.SpikeIdMat(:,[1,2]), [Fig.Current.DateNo, Fig.Current.ChannelNo],'rows')),3))'),'%s');
            set(Fig.Pannel.InputsH(4),'string', Data.Cells,'value',Fig.Current.CellNo);
           	PlotCellData;
            
        case 4
            Fig.Current.CellNo = get(hObj,'value');
            PlotCellData;
    end
end
