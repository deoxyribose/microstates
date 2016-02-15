function [OUTEEG, com] = plotsegmentation(INEEG,frame_start,frame_end);

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

if ~isstruct(INEEG) % in this case, assume INEEG is a T x n matrix of segmentations
    OUTEEG.data = INEEG;
    OUTEEG.pnts = size(INEEG,1);
    OUTEEG.K = numel(unique(INEEG(:,1)));
    OUTEEG.srate = OUTEEG.pnts;
else
    OUTEEG = INEEG;
end
if ~isfield(OUTEEG, 'K')
    error('You must estimate the microstates first');
end
if ~isfield(OUTEEG, 'gtm')
    OUTEEG.gtm = double(std(OUTEEG.data,[],1));
    OUTEEG.gtmtype = 1;
end


if OUTEEG.pnts > 1000 & nargin == 1
	promptstr    = { 'Plot from frame nr:', 'Plot to frame nr:' };
	inistr       = { '1',num2str(min(OUTEEG.srate,OUTEEG.pnts))};
	result       = inputdlg( promptstr, 'Plot topographic measures', 1,  inistr);
	if length( result ) == 0 return; end;

	frame_start   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
    frame_end   	 = eval( [ '[' result{2} ']' ] ); % the brackets allow to process matlab arrays
    if frame_end > OUTEEG.pnts
        frame_end = OUTEEG.pnts;
    end
elseif nargin == 1
    frame_start = 1;
    frame_end = numel(OUTEEG.gtm);
end;

% a plot of the gfp, where areas between peaks are colored according to
% OUTEEG.Z
figure('name','Microstate segementation','Position', [100, 100, 900, 100]);
x = OUTEEG.times;
if OUTEEG.gtmtype == 1
    y = OUTEEG.gtm;
else
    y = -1*OUTEEG.gtm';
end
c = OUTEEG.Z;
lengths = [size(x,2) size(y,2) size(c,2)];
if range(lengths) ~= 0 % check for equal lengths
    min_length = min(lengths);
    x = x(1:min_length);
    y = y(1:min_length);
    c = c(1:min_length);
    frame_end = min(min_length,frame_end);
end
x = x(frame_start:frame_end);
y = y(frame_start:frame_end);
if size(c,1) > 1
    for i=1:OUTEEG.K
        subplot(OUTEEG.K,1,i)
        c_i = c(i,frame_start:frame_end)+1;
        surface([x;x],[zeros(size(x));y],[zeros(size(x));zeros(size(x))],[c_i;c_i],'facecol','no','edgecol','flat','linew',2);
        grid on;
        xlabel('ms')
        ylabel('GFP')
        %colorbar('Ticks',1:OUTEEG.K)
    end
else
    c = c(frame_start:frame_end);
    surface([x;x],[zeros(size(x));y],[zeros(size(x));zeros(size(x))],[c;c],'facecol','no','edgecol','flat','linew',2);
    grid on;
    xlabel('ms')
    ylabel('GFP')
    colorbar('Ticks',1:OUTEEG.K)
    %legend(num2str(1:OUTEEG.K),num2str(1:OUTEEG.K))
end
com = sprintf('pop_plotsegmentation( %s, %d, [%s] );', inputname(1), frame_start, frame_end);

