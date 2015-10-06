function show_clusters(signal,clustering,meanGMD,K)
if nargin < 2 % then signal = EEG
    EEG = signal;
    if isfield(EEG, 'idx')
        clustering = EEG.idx;
        K = EEG.nMCRSTS;
    else
        clustering = ones(1,size(EEG.data,2));
    end
    signal = EEG.data;
end
if size(signal,2) > 5000
    signal = signal(:,1:5000)';
    if ~isfield(EEG, 'idx')
        clustering = ones(1,5000);
    else
        clustering = EEG.idx(1:5000);
    end
else
    signal = signal';
end
if nargin < 5
    % standardize -> pca -> plot 3 biggest var pc's
    % zero mean and scale
    X=bsxfun(@rdivide,bsxfun(@minus,signal,mean(signal,2)),std(signal));
    % visualize clustering by PCA
    [U,S,~]=svdecon(X);
    signal = U*S;
    vars = cumsum(diag(S).^2)/trace(S.^2);
end


%K0 = find(vars<0.99);
%K0 = K0(end);
figure('name', '3d plot of clustering','Position',[100,100,900,900])
%for k=1:K0,
scatter3(signal(:,1),signal(:,2),signal(:,3),100*(signal(:,4)+abs(min(signal(:,4)))+0.1),clustering,'.')
if nargin == 4
    title(['variance explained = ',num2str(vars(3)), ', mean GMD = ',num2str(meanGMD),'   |   k = ', num2str(K)])
else
    title(['variance explained = ',num2str(vars(3))])
end
%hold on;
%text(V(:,1),V(:,3),int2str(k));
%end
end