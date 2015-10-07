function [OUTEEG, com] = global_topographic_measures( INEEG, typeproc, frame_start, frame_end );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
OUTEEG = INEEG;
% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help global_topographic_measures;
	return;
end;	

% pop up window
% -------------
if nargin < 3
	promptstr    = { 'Plot from frame nr:', 'Plot to frame nr:' };
	inistr       = { '1',num2str(OUTEEG.srate) };
	result       = inputdlg( promptstr, 'Plot topographic measures', 1,  inistr);
	if length( result ) == 0 return; end;

	frame_start   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
    frame_end   	 = eval( [ '[' result{2} ']' ] ); % the brackets allow to process matlab arrays
    if frame_end > OUTEEG.pnts
        frame_end = OUTEEG.pnts
    end
end;

% call function sample either on raw data or ICA data
% ---------------------------------------------------
if typeproc == 1
	dat = OUTEEG.data(:,frame_start:frame_end);
    OUTEEG.gtm = std(dat,[],1);
    %tend = 500;
    figure;
    x = OUTEEG.times(frame_start:frame_end);
    x = x(2:end);
    gfp = OUTEEG.gtm(2:end);
    gmd = GMD(dat(:,1:end-1),dat(:,2:end),OUTEEG.nbchan);
    R = corrcoef(gfp,gmd);
    subplot(211)
    suptitle(['Global Topographic Measures ', num2str(R(1,2))])
    plot(x, gfp)
    title('Global Field Power')
    subplot(212)
    plot(x, gmd);
    xlabel('ms');
    %ylabel('Standard deviation across channes');
    title('Global Map Dissimilarity')
    
else
	if ~isempty( OUTEEG.icadata )
		sample( OUTEEG.icadata );
	else
		error('You must run ICA first');
	end;	
end;
% return the string command
% -------------------------
com = sprintf('global_topographic_measures( %s, %d, [%s] );', inputname(1), typeproc, int2str(frame_start), int2str(frame_end));

return;


end

