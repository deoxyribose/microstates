function [OUTEEG, com] = pop_get_clean_topographies( INEEG, typeproc);

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
OUTEEG = INEEG;
% display help if not enough arguments
% ------------------------------------

% pop up window
% -------------
if nargin < 2
	promptstr    = { 'Choose topographic measure' };
	inistr       = { '1 for gfp, 2 for gmd' };
	result       = inputdlg( promptstr, 'Title of window', 1,  inistr);
	if length( result ) == 0 return; end;

	OUTEEG.gtmtype	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
else 
    OUTEEG.gtmtype = typeproc;
end;
switch OUTEEG.gtmtype;
    case 1  % Global Field Power
        OUTEEG.gtm = double(std(OUTEEG.data,[],1));
        OUTEEG.gtm = OUTEEG.gtm(2:end);
        [~, OUTEEG.peakidx] = findpeaks(OUTEEG.gtm,'MinPeakDistance',OUTEEG.srate/20); %min peak dist 5 ms
        figure('name','GFP peaks');plotSignalWithPeaks(OUTEEG.times(1:OUTEEG.srate)/2,OUTEEG.gtm(1:OUTEEG.srate),OUTEEG.peakidx(OUTEEG.peakidx<OUTEEG.srate));
        figure('name','GFP peaks');plotSignalWithPeaks(OUTEEG.times(1:OUTEEG.srate*10)/2,OUTEEG.gtm(1:OUTEEG.srate*10),OUTEEG.peakidx(OUTEEG.peakidx<OUTEEG.srate*10));
    case 2
        OUTEEG.gtm = -1*GMD(OUTEEG.data(:,1:end-1),OUTEEG.data(:,2:end),OUTEEG.nbchan);
        %tmp = movingAverage(OUTEEG.gtm,5)./movingAverage(OUTEEG.gtm,100);
        %subplot(211)
        %plot(double(OUTEEG.gtm(1:30000)))
        %subplot(212)
        %plot(tmp(1:30000))
        %[~, OUTEEG.peakidx] = findpeaks(double(movingAverage(OUTEEG.gtm(1:10000),20)),'MinPeakHeight',0.02,'MinPeakDistance',30);
        [~, OUTEEG.peakidx] = findpeaks(double(OUTEEG.gtm),'MinPeakDistance',OUTEEG.srate/100); %'MinPeakHeight',-0.1,
        figure;plotSignalWithPeaks(OUTEEG.times(1:OUTEEG.srate*10)/2,OUTEEG.gtm(1:OUTEEG.srate*10),OUTEEG.peakidx(OUTEEG.peakidx<OUTEEG.srate*10));
    otherwise
        subplot(311)
        plot(OUTEEG.gfp(1:5000))
        subplot(312)
        plot(movingAverage(diff(OUTEEG.gfp(1:5000)),20))
        subplot(313)
        plot(movingAverage(OUTEEG.gmd(1:5000),20))
end
disp(['Found ', num2str(length(OUTEEG.peakidx)), ' peaks, which is ', num2str(length(OUTEEG.peakidx)/OUTEEG.pnts*100), ' percent of all data']);
return;