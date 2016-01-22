function [mse] = plot_true_vs_estimate(w,x,s,W,X,Ms,chan_idx,time_idx,chanlocs)
%w = w(chan_idx,:);
Wtmp = w;
%Wtmp(chan_idx,:) = W;
x = x(:,time_idx);
s = s(:,time_idx);
[K,T] = size(x);
figure('name','Simulated microstate spatial topography maps');
tmp = corr(W,w(chan_idx,:));
[~,I] = max(abs(tmp),[],2);
polarity = round(diag(tmp(:,I)));
for i=1:K
    subplot(2,K,i)
    topoplot(w(:,i),chanlocs,'plotchans',chan_idx,'electrodes','off','style','map');
    subplot(2,K,i+K)
    topoplot(polarity(i)*Wtmp(:,i),chanlocs,'plotchans',chan_idx,'electrodes','off','style','map');
end
figure
subplot(4,1,1)
plot(repmat(1:T,K,1)',(s.*x)')
title('True s.*x')
subplot(4,1,2)
title('Estimated Ms.*X')
plot(repmat(1:T,K,1)',(Ms.*X)')
subplot(4,1,3)
imagesc(s)
title('True s')
subplot(4,1,4)
title('Estimated S')
imagesc(Ms)
end