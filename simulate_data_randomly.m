function [w,x,s,yrec,p0] = simulate_data_randomly(J,K,T,draw,chanlocs)
% forward model
% load('/home/frans/documents/Microstates/EEG_faces/SPMgainmatrix_aceMdspm8_faces_run1_1_c.mat')
w=rand(J,K);
w=bsxfun(@rdivide,bsxfun(@minus,w,mean(w,1)),std(w,[],1)); % normalize
if draw
    chan_idx = 1:(128/J):128;
    dims = factor(K);
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
    figure;
    for i=1:K
        subplot(nrow,ncol,i)
        topoplot(w(:,i),chanlocs(chan_idx),'electrodes','off','style','map');
    end
end
x = rand(K,T);

s = zeros(K,T);
for k=1:K
    s(k,:) = abs(sin(0.005*(1:T)))>(k-1)/K;
    s(k,abs(sin(0.005*(1:T)))>k/K) = 0;
end
if draw
    subplot(2,1,1)
    imagesc(s)
    subplot(2,1,2)
    plot(repmat(1:T,K,1)',(x.*s)')
end
[~,I] = max(s);
p0 = 1-sum(abs(diff(I)))/T;
yrec = w*(x.*s);
end