function [w,x,s,yrec,p0] = simulate_data(J,K,T,G,draw,chanlocs)
% forward model
% load('/home/frans/documents/Microstates/EEG_faces/SPMgainmatrix_aceMdspm8_faces_run1_1_c.mat')
chan_idx = 1:(128/J):128;
G = G(chan_idx,:);
w = rand(size(G,2),K);
w(w<.95) = 0;
w = G*w;                  % transform to scalp "measurements"
% w = rand(J,K);
w=bsxfun(@rdivide,bsxfun(@minus,w,mean(w,1)),std(w,[],1)); % normalize
if draw
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
% random frequency, random phase, random amplitude
%x = Mu + rand(K,T)*alpha;

freqs1 = (1+rand(1,1).^([1:K])')/10;
freqs2 = (1+rand(1,1).^([1:K])')/10;
phase1 = rand(K,1);
phase2 = rand(K,1);
amp1 = [1:K]'*50;
amp2 = [1:K]'*50;
offset = repmat(rand(K,1)*100,1,T);
x = 0.1*bsxfun(@times,sin(bsxfun(@plus,bsxfun(@times,repmat(1:T,K,1),freqs1),phase2)),amp1) ...
    + 0.1*bsxfun(@times,cos(bsxfun(@plus,bsxfun(@times,repmat(1:T,K,1),freqs2),phase2)),amp2) + offset;


s = zeros(K,T);
for k=1:K
    s(k,:) = abs(sin(0.002*(1:T)))>(k-1)/K;
    s(k,abs(sin(0.002*(1:T)))>k/K) = 0;
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