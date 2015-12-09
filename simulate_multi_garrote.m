J = 128;
K = 3;
T = 512;
alpha = 0.8;
beta = 6;

trans = rand(K,K)+diag((1:K));
trans = bsxfun(@rdivide,trans,sum(trans,2));

% forward model
load('/home/frans/documents/Microstates/EEG_faces/SPMgainmatrix_aceMdspm8_faces_run1_1_c.mat')

%w = randi(size(G,2),1,K); % choose source randomly
w = [500, 5501, 7502];
w = arr2mat(w,size(G,2)); % make it a matrix
w = G*w;                  % transform to scalp "measurements"
w=bsxfun(@rdivide,bsxfun(@minus,w,mean(w,2)),std(w,[],2)); % normalize
figure('name','Simulated microstate spatial topography maps');
for i=1:K
    subplot(5,3,i)
    topoplot(w(:,i),EEG.chanlocs,'electrodes','off','style','map');
end
% random frequency, random phase, random amplitude
Mu = 0.1*bsxfun(@times,sin(bsxfun(@plus,bsxfun(@times,repmat(1:T,K,1),randi(K,1,K)'/10^2),randi(20,K,1))) ...
    + cos(bsxfun(@plus,bsxfun(@times,repmat(1:T,K,1),randi(K,1,K)'/10^2),randi(20,K,1))),randi(20,K,1));
x = Mu + rand(K,T)*alpha;
figure
plot(repmat(1:T,K,1)',x')
s = zeros(1,T);
s(1) = 1;
for i = 2:T
    s(i) = datasample(1:K,1,'Weights',trans(s(i-1),:));
end
    
s = arr2mat(s);
figure
imagesc(s)
figure
plot(repmat(1:T,K,1)',(x.*s)')


Y = w*(x.*s) + rand(J,T)*beta^(-1/2);
figure
plot(repmat(1:T,J,1)',Y(1:J,:)')
figure('name','Simulated microstate spatial topography maps');
for i=1:K
    subplot(5,3,i)
    topoplot(w(:,i),EEG.chanlocs,'electrodes','off','style','map');
end

%%
%W = bsxfun(@rdivide,Y*(Mu.*M)',sum((Mu.^2+sigma2).*M,2)');
W = bsxfun(@rdivide,Y*(Mu.*s)',sum((Mu.^2+sigma2).*s,2)');
figure
plot(w,W,'*')
figure
plot(w,bsxfun(@rdivide,bsxfun(@minus,W,mean(W,1)),std(W,[],1)),'*')
axis equal

WcircW = w.*w;
diagWW=sum(WcircW,1);
WY = w'*Y;

sigma2 = 1./((1/alpha^2)+beta*bsxfun(@times,s,diagWW')); % what is sigma2 supposed to be?

Mu_hat = beta*(WY.*s).*sigma2;
tmp = Mu.*s;
figure
plot(tmp(:,:)',Mu_hat(:,:)','*')

beta_hat = (J*T)/sum(sum(Y.*Y- 2*Y.*(w*(Mu.*s)) +WcircW*((Mu.^2 + sigma2).*s)));

M=softmax(1*beta*( (WY).*Mu - 0.5*diag(diagWW)*(Mu.^2 + sigma2)));
[~,Ms]=max(M,[],1);
Ms=arr2mat(Ms,K);
figure
imagesc(Ms)
%% Example
EEG = pop_loadset('filename','simulated3.set','filepath','/home/frans/MATLAB packages/eeglab13_4_4b/plugins/microstates0.1/');
%load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat')
%EEG.chanlocs = chanlocs(1:30);
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5); % auto-clean plugin
EEG=pop_chanedit(EEG, 'load',{'/home/frans/documents/Microstates/bdf-2-M-Neuch/biosemi128.xyz' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
eeglab('redraw')
%%
EEG = pop_loadset('filename','simulated_3_msts.set','filepath','/home/frans/MATLAB packages/eeglab13_4_4b/plugins/microstates0.1/');
%load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat')
%EEG.chanlocs = chanlocs(1:30);
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5); % auto-clean plugin
EEG=pop_chanedit(EEG, 'load',{'/home/frans/documents/Microstates/bdf-2-M-Neuch/biosemi128.xyz' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
eeglab('redraw')
%%
EEG = pop_loadset('filename','simulated.set','filepath','/home/frans/MATLAB packages/eeglab13_4_4b/plugins/microstates0.1/');
%load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat')
%EEG.chanlocs = chanlocs(1:30);
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5); % auto-clean plugin
EEG=pop_chanedit(EEG, 'load',{'/home/frans/documents/Microstates/bdf-2-M-Neuch/biosemi128.xyz' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
eeglab('redraw')
