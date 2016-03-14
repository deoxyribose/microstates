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

[J,T,nepochs] = size(Y);
c = 1;
for n = [1,4,10,21,43,86]
    % divide Y into cells of JxTxn, and one cell with the number of epochs
    % remaining after division by n, and mean over the epochs that have
    % equal number of events
    this_isnt_even_my_final_form = mat2cell(Y, [J],[T],[ones(1,floor(nepochs/n))*n rem(nepochs,n)]);
    spm_erps{c} = cell2mat(cellfun(@(x)mean(x,3),this_isnt_even_my_final_form(1:end-1),'UniformOutput',false));
    c = c+1;
end
nerps = size(spm_erps,3);

nsnrs = 6;
nexperiment = 50;
clustering_algorithms = [4,6];
nalgo = numel(clustering_algorithms);
allidxs = [3, nalgo, nsnrs,]
for nmicrostates = [2,3,4]
    for c=1:nalgo
        for j = 1:nsnrs
            tmp = spm_erps{j};
            for k = 1:size(spm_erps{j},3)
                EEG.data = tmp(:,:,k);
                parfor i=1:nexperiment
                    OUTEEG = pop_getmicrostates(EEG,'Kfrom',nmicrostates,'Kto',nmicrostates,'clustering_algorithm',clustering_algorithms(c),'chronos',0,'MULTI',60,'verbose',0);
                    labelings(:,c,j,k,i) = OUTEEG.Z;
                end
            end
        end
    end
end