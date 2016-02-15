%% Load data
warning('off','all')
addpath(genpath('/home/frans/MATLAB packages/'));
rmpath(genpath('/home/frans/MATLAB packages/eeglab13_5_4b/plugins/Biosig2.88/'));
cd('/home/frans/MATLAB packages/eeglab13_5_4b/plugins/microstates0.1')
%load('newerps.mat')
%STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[]; % clear any previous study
% Original-my17ss-1-concat-650ms.eph
EEG = pop_loadset('filename','violaine.set','filepath','/home/frans/documents/Microstates/');
%load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat')
%EEG.chanlocs = chanlocs(1:30);
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5); % auto-clean plugin
EEG=pop_chanedit(EEG, 'load',{'/home/frans/documents/Microstates/bdf-2-M-Neuch/biosemi128.xyz' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
eeglab('redraw')

load bigfif.mat
%% Model selection using MI and entropy?
% size(MIs)
% 
% ans =
% 
%    100     3     2   100
%   nbeta   nk   algo  inits
% size(Hs)
% 
% ans =
% 
%    100     3     2   100
[~,Z_true] = max(s,[],1);
entropy(Z_true) % = mutualinfo(Z_true,Z_true), overlapping information sets
plot(beta,MIs(:,:,1,1),beta,Hs(:,:,1,1))

figure; 
subplot(3,1,1);
plot(beta,MIs(:,:,1,1));
legend('4','7','10');
subplot(3,1,2);
plot(beta,Hs(:,:,1,1));
legend('4','7','10');
subplot(3,1,3);
plot(beta,1.5*MIs(:,:,1,1)-(Hs(:,:,1,1))); 
legend('4','7','10');

figure;
c=1;
for i=1:3
    for j=1:2
        subplot(3,4,c)
        imagesc(squeeze(MIs(:,i,j,:)))
        ylabel('SNR')
        xlabel('Model #')
        title(['I(Z_true,Z_hat),K=',num2str(1+i*3)])
        axis('equal')
        axis('tight')
        caxis manual
        caxis([0 3]);
        colorbar()
        c = c+1;
        subplot(3,4,c)
        imagesc(squeeze(Hs(:,i,j,:)))
        ylabel('SNR')
        xlabel('Model #')
        title(['H(Z_est),K=',num2str(1+i*3)])
        axis('equal')
        axis('tight')
        caxis manual
        caxis([0 3]);
        colorbar()
        c = c+1;
    end
end


figure;
c=1;
for i=1:3
    for j=1:2
        subplot(3,4,c)
        [bandwidth,density,xmesh,cdf]=kde(MIs(:,:,1,:));
        plot(xmesh,density)
        ylabel('P')
        xlabel('MI')
        title(['I(Z_true,Z_hat),K=',num2str(1+i*3)])
%         axis('equal')
        axis('tight')
        c = c+1;
        subplot(3,4,c)
        [bandwidth,density,xmesh,cdf]=kde(MIs(:,:,2,:));
        plot(xmesh,density)
        ylabel('P')
        xlabel('MI')
        title(['H(Z_est),K=',num2str(1+i*3)])
%         axis('equal')
        axis('tight')
        c = c+1;
    end
end



% argmin softmax(c1*MIs+c2*Hs)