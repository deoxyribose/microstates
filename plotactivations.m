function [OUTEEG, com] = plotmicrostates( INEEG);

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
OUTEEG = INEEG;
% display help if not enough arguments
% ------------------------------------
% Check whether getmicrostates has been run succesfully.
if ~isfield(OUTEEG, 'K') 
    error('You must estimate the microstates first');
end
if min(size(OUTEEG.Z)) == 1
    Z = arr2mat(OUTEEG.Z,OUTEEG.K);
    usedmsts = unique(OUTEEG.Z);
    OUTEEG.K = length(usedmsts);
else
    Z = OUTEEG.Z;
end
figure
plot(repmat(1:OUTEEG.pnts,OUTEEG.K,1)',(Z.*OUTEEG.A)')
legend('show')

end