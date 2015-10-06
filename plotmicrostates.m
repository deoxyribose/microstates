function [OUTEEG, com] = plotmicrostates( INEEG, typeproc, param3 );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
OUTEEG = INEEG;
% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_plotmicrostates;
	return;
end;

% check for channel locations
if ~isfield(OUTEEG, 'chanlocs') 
    error('Load channel locations first');
end

% Check whether getmicrostates has been run succesfully.
if ~isfield(OUTEEG, 'nMCRSTS') 
    error('You must estimate the microstates first');
end

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

figure('name','Microstate spatial topography maps');
for i=1:OUTEEG.nMCRSTS
    subplot(nrow,ncol,i)
    %lo = min(min(OUTEEG.A));
    %hi = max(max(OUTEEG.A));
    %topoplot(OUTEEG.A(:,i),OUTEEG.chanlocs,'electrodes','off','style','map','maplimits',[lo,hi]);%,'plotchans',chanlocs(1:30).labels);
    topoplot(OUTEEG.A(:,i),OUTEEG.chanlocs,'electrodes','off','style','map');%,'plotchans',chanlocs(1:30).labels);
    if OUTEEG.clustering_algorithm == 1 % if ICA has been used
        title(strvcat(['Microstate ',num2str((OUTEEG.nMCRSTS+1)-i)],num2str(OUTEEG.vars(i))))
    else
        title(['Microstate ',num2str(i)])
    end
end
btn = uicontrol('Style', 'pushbutton', 'String', 'Close','Units','normalized',...
        'Position', [0.85 0.01 0.05 0.05],...
        'Callback', 'delete(gcf)');     
com = sprintf('plotmicrostates( %s, %d, [%s] );', inputname(1), typeproc, int2str(param3));

end