function microstate_statistics( INEEG, display_in_command_window);

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
% if the user press the cancel button
OUTEEG = INEEG;
% display help if not enough arguments
% ------------------------------------
if nargin < 2
    promptstr    = { 'Display in MATLAB command window? (yes, no) = (1,0)' };
    inistr       = { '0' };
    result       = inputdlg( promptstr, 'Title of window', 1,  inistr);
    if length( result ) == 0 return; end;
    display_in_command_window   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
end;

if ~isfield(OUTEEG, 'nMCRSTS')
    error('You must estimate the microstates first');
end
% Determine the kind of analysis from paramters from get_microstates
switch OUTEEG.subset
    case 1
        subset = ['GFP-peak topographies, corresponding to ', num2str(length(OUTEEG.peakidx)/OUTEEG.pnts)];
    case 2
        subset = 'all';
    case 3
        subset = num2str(length(OUTEEG.peakidx)/OUTEEG.pnts*100);
end
switch OUTEEG.clustering_algorithm
    case 1
        clustering_algorithm = 'ICA';
    case 2
        clustering_algorithm = 'kmeans';
    case 3
        clustering_algorithm = 'agglomerative clustering';
end

analysis_description = ['Estimated ', num2str(OUTEEG.nMCRSTS), ' microstates, using ', clustering_algorithm, ' on ', subset, ' of the data'];
        
A = OUTEEG.idx';
J=find(diff([A(1)-1; A]));
microstate_durations=[A(J), diff([J; numel(A)+1])];

avg_nframes = zeros(OUTEEG.nMCRSTS,1);
nsegments = zeros(OUTEEG.nMCRSTS,1);
all_segment_durations = cell(20,1);
for microstate = 1:OUTEEG.nMCRSTS
    durations_of_current_microstate = microstate_durations(microstate_durations(:,1)==microstate,2);
    all_segment_durations{microstate} = [all_segment_durations{microstate} durations_of_current_microstate];
    avg_nframes(microstate) = mean(durations_of_current_microstate);
    nsegments(microstate) = numel(durations_of_current_microstate);
end
totalnsegments = sum(nsegments);
total_duration_in_frames = round(avg_nframes.*nsegments,0);
% Calculate number of rows and cols in the subplot figure
dims = factor(OUTEEG.nMCRSTS);
if size(dims,2)==4
    nrow = dims(1)*dims(2);
    ncol = dims(3)*dims(4);    
elseif size(dims,2)==3
    nrow = dims(1)*dims(2);
    ncol = dims(3);
elseif size(dims,2)==2
    nrow = dims(1);
    ncol = dims(2);
else
    nrow = 1;
    ncol = dims(1);
end
% histograms of durations
figure(1);
hist(cell2mat(all_segment_durations),100)
title(['log-log histogram of microstate durations, mean is ',num2str(mean(avg_nframes))])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
btn = uicontrol('Style', 'pushbutton', 'String', 'Close','Units','normalized',...
        'Position', [0.85 0.01 0.05 0.05],...
        'Callback', 'delete(gcf)');     
% histograms of durations per microstate
% promptstr    = { 'Histograms or boxplots? (1, 0)' };
% inistr       = { '0' };
% result       = inputdlg( promptstr, 'Title of window', 1,  inistr);
% if length( result ) == 0 return; end;
% histograms = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
histograms = 1;
boxplots = 0;
if histograms
    figure(2);
    c = 1;
    for i = 1:nrow
        for j = 1:ncol
            subplot(nrow,ncol,c)
            hist(all_segment_durations{c},80)
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
            %xlim([-5,5])
            c = c+1;
        end
    end
else
    figure(3);
    c = 1;
    for i = 1:nrow
        for j = 1:ncol
            subplot(nrow,ncol,c)
            boxplot(all_segment_durations{c})
            c = c+1;
        end
    end
end
btn = uicontrol('Style', 'pushbutton', 'String', 'Close','Units','normalized',...
        'Position', [0.85 0.01 0.05 0.05],...
        'Callback', 'delete(gcf)');     
f = figure(4);
%title(analysis_description)
cnames = {'Microstate nr','mean duration (ms)','nr. of segments','total duration (s)','percentage of signal'};
format longG
info = [(1:OUTEEG.nMCRSTS)' round(avg_nframes*(1000/OUTEEG.srate),0) nsegments round(total_duration_in_frames/OUTEEG.srate) round(total_duration_in_frames/OUTEEG.pnts*100,1)];
transmat = gettransitionMatrix(OUTEEG);

if ~display_in_command_window
    t = uitable(f,'Data',info,'ColumnName',cnames);
    t.Position(1) = t.Extent(1);
    t.Position(2) = t.Extent(2);
    t.Position(3) = t.Extent(3);
    t.Position(4) = t.Extent(4);
    % from http://undocumentedmatlab.com/blog/uitable-sorting
    jscrollpane = findjobj(t);
    jtable = jscrollpane.getViewport.getView;
    % Now turn the JIDE sorting on
    jtable.setSortable(true);		% or: set(jtable,'Sortable','on');
    jtable.setAutoResort(true);
    jtable.setMultiColumnSortable(true);
    jtable.setPreserveSelectionsAfterSorting(true);
    btn = uicontrol('Style', 'pushbutton', 'String', 'Close','Units','normalized',...
        'Position', [0.85 0.01 0.05 0.05],...
        'Callback', 'delete(gcf)');     
else
    if numel(microstate_durations)<50
        microstate_durations
    end
    disp(cnames);
    disp(info);
    disp('Correlation matrix:');    
    disp(corrcoef(OUTEEG.A));
    disp('Transition matrix:');
    disp(transmat);
end
    com = sprintf('microstate_statistics( %s, %d, [%s] );', inputname(1));
end