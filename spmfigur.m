warning('off','all')
addpath(genpath('/home/frans/MATLAB packages/'));
rmpath(genpath('/home/frans/MATLAB packages/eeglab13_5_4b/plugins/Biosig2.88/'));
cd('/home/frans/MATLAB packages/eeglab13_5_4b/plugins/microstates0.1')
load('EEGfaces_JxTxn_epochs.mat') % size(Y) = 128 * 1639 * 172

STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[]; % clear any previous study
EEG = pop_loadset('filename','EEGfaces_epoch_2_Vmicrostates.set','filepath','/home/frans/MATLAB packages/eeglab13_5_4b/plugins/microstates0.1');
EEG=pop_chanedit(EEG, 'load',{'/home/frans/Dropbox/Microstates/chanlocs.sfp' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
eeglab('redraw')
condition = textread('/home/frans/Dropbox/Microstates/face_scrambled_face_SPM_data/condition_labels.txt','%s');
Yall = Y;
Y = Yall(:,:,strcmp(condition,'faces'));
Y = faces_baseline_correction(Y);
[J,T,nepochs] = size(Y);
c = 1;
for n = [5,10,21,43]
    % divide Y into cells of JxTxn, and one cell with the number of epochs
    % remaining after division by n, and mean over the epochs that have
    % equal number of events
    cell_of_erps = mat2cell(Y, [J],[T],[ones(1,floor(nepochs/n))*n rem(nepochs,n)]);
    spm_erps{c} = cell2mat(cellfun(@(x)mean(x,3),cell_of_erps(1:end-1),'UniformOutput',false));
    c = c+1;
end
nnmicrostates = 1;
nsnrs = 4; % 4
nexperiment = 20; % 20
clustering_algorithms = [4,6];
nalgo = numel(clustering_algorithms);
all_erps = cat(3,spm_erps{:}); % [17,8,4,2]
nerps = size(all_erps,3);
EEG.data = all_erps;

allidx = [nnmicrostates,nalgo,nerps,nexperiment];
parfor idx = 1:prod(allidx)
    [nmicrostates,c,j,i] = ind2sub(allidx,idx);
    OUTEEG = pop_getmicrostates(EEG,'Kfrom',nmicrostates+1,'Kto',nmicrostates+1,'clustering_algorithm',clustering_algorithms(c),'chronos',0,'MULTI',2,'verbose',0,'j',j); %MULTI = 60
    labelings(:,idx) = OUTEEG.Z;
end
save('labelings5exps.mat','labelings')
% for nmicrostates = 2
%     for c=1:nalgo
%         for j = 1:nerps
%             EEG.data = all_erps(:,:,j);
%             parfor i=1:nexperiment
%                 OUTEEG = pop_getmicrostates(EEG,'Kfrom',nmicrostates,'Kto',nmicrostates,'clustering_algorithm',clustering_algorithms(c),'chronos',0,'MULTI',60,'verbose',0);
%                 labelings(:,c,j,i) = OUTEEG.Z;
%             end
%         end
%     end
% end


% for nmicrostates = 2
%     for c=1:nalgo
%         for j = 1:nsnrs
%             tmp = spm_erps{j};
%             for k = 1:size(spm_erps{j},3)
%                 EEG.data = tmp(:,:,k);
%                 parfor i=1:nexperiment
%                     OUTEEG = pop_getmicrostates(EEG,'Kfrom',nmicrostates,'Kto',nmicrostates,'clustering_algorithm',clustering_algorithms(c),'chronos',0,'MULTI',60,'verbose',0);
%                     labelings(:,c,j,k,i) = OUTEEG.Z;
%                 end
%             end
%         end
%     end
% end