function [OUTEEG, com] = merge_correlated_microstates(INEEG);

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
OUTEEG = INEEG;
corrmat = abs(corr(OUTEEG.W));
meancorr = sort(mean(corrmat-diag(diag(corrmat)),1));
figure;subplot(2,1,1);plot(meancorr);title('Mean correlation of a microstate to the other microstates \n Click on preferred threshold correlation. Hit Enter.');axis('tight');subplot(2,1,2);plot(diff(meancorr));axis('tight');title(' Everything correlated more than the threshold will be merged');
[x,y] = ginput;
x = x(end);
y = y(end);
close
% display help if not enough arguments
% ------------------------------------
[~, above_threshold_idx] = find(meancorr>y);
idx = [1:size(OUTEEG.W,2)];
tmp = setdiff(idx,above_threshold_idx);
kmax = max(tmp);
if ~all(1:kmax == sort(tmp))
    keyboard
end
kmax = kmax+1;
OUTEEG.K = kmax;
OUTEEG.W(:,kmax) = mean(OUTEEG.W(:,above_threshold_idx),2);
OUTEEG.W = OUTEEG.W(:,1:kmax);
OUTEEG.A(kmax,:) = mean(OUTEEG.A(above_threshold_idx,:),1);
OUTEEG.A = OUTEEG.A(1:kmax,:);
OUTEEG.Z(OUTEEG.Z>kmax-1) = kmax;
end