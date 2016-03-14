%%
warning('off','all')
addpath(genpath('/home/frans/MATLAB packages/'));
rmpath(genpath('/home/frans/MATLAB packages/eeglab13_5_4b/plugins/Biosig2.88/'));
cd('/home/frans/MATLAB packages/eeglab13_5_4b/plugins/microstates0.1')
load('EEGfaces_JxTxn_epochs.mat') % size(Y) = 128 * 1639 * 172
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[]; % clear any previous study
EEG = pop_loadset('filename','EEGfaces_epoch_1_polymicrostates.set','filepath','/home/frans/MATLAB packages/eeglab13_5_4b/plugins/microstates0.1');
EEG=pop_chanedit(EEG, 'load',{'/home/frans/Dropbox/Microstates/chanlocs.sfp' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
eeglab('redraw')
%% make ERPs from n epochs
% two ways to do this: either move in non-overlapping windows of n, which is standard, or
% overlapping windows of n, which gives "more" data
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
%% Grand ERP
Y = mean(Y,3);
Y=bsxfun(@rdivide,bsxfun(@minus,Y,mean(Y,2)),std(Y,[],2));

imagesc(corr(Y));colorbar()
[U,S,~]=svdecon(Y');
signal = U*S;
figure;scatter(signal(:,1),signal(:,2),15,[1:size(Y,2)],'.');colorbar()

% permute channels so that nearby idx's are correlated
largeabscorr = abs(corr(Y'))>.5;
figure;
spy(largeabscorr)
chan_idx = amd(largeabscorr);
figure;
spy(abs(corr(Y(chan_idx,1:500)'))>.5)
Y = Y(chan_idx,:);
%% Performance vs SNR
%% Majority vote segmentation of grand ERP in 50 experiments
nsnrs = 6;
nexperiment = 50;
clustering_algorithms = [4,6];
nalgo = numel(clustering_algorithms);
%recon_error = zeros(nalgo,nexperiment);
%algoname = {'k_means','agglomerative','n_microstates','variational_microstates','polymicrostates'};
%OUTEEGs = cell(nalgo);
for nmicrostates = 2%[3,7]
    for c=1:nalgo
        for j = 4:4%nsnrs
            tmp = spm_erps{j};
            for k = 1:size(spm_erps{j},3)
                EEG.data = tmp(:,:,k);
                parfor i=1:nexperiment
                    %std(EEG.data(:,1),[],1)
                    % 2 = kmeans, 3 = agglo, 4 = n-microstates, 5 = dpmeans
                    % 6 = variational microstates, 7 = polymicrostates
                    OUTEEG = pop_getmicrostates(EEG,'Kfrom',nmicrostates,'Kto',nmicrostates,'clustering_algorithm',clustering_algorithms(c),'chronos',0,'MULTI',60,'verbose',0);
                    %OUTEEGs{c} = OUTEEG;
                    labelings(:,c,j,k,i) = OUTEEG.Z;
                end
            end
        end
    end
end

save('labelings5.mat','labelings');
max_snr_Z_n = squeeze(labelings(:,1,6,1,:));
max_snr_Z_v = squeeze(labelings(:,2,6,1,:));

max_snr_Z_n = squeeze(labelings(:,1,6,2,:));
max_snr_Z_v = squeeze(labelings(:,2,6,2,:));

figure;
subplot(121)
imagesc(max_snr_Z_n')
subplot(122)
imagesc(max_snr_Z_v')

[est_perm,signs]=soft_align(arr2mat(max_snr_Z_n(:,end)',nmicrostates),arr2mat(max_snr_Z_v(:,1)',nmicrostates));
%majority_vote_n = 
%% Test reconstruction of different algorithms
nexperiment = 1;
clustering_algorithms = [4,6,7];
nalgo = numel(clustering_algorithms);
recon_error = zeros(nalgo,nexperiment);
EEG.data = Y;
algoname = {'k_means','agglomerative','n_microstates','variational_microstates','polymicrostates'};
OUTEEGs = cell(nalgo);
for nmicrostates = 2%[3,7]
    for c=1:nalgo
        for i=1:nexperiment
            %std(EEG.data(:,1),[],1)
            % 2 = kmeans, 3 = agglo, 4 = n-microstates, 5 = dpmeans
            % 6 = variational microstates, 7 = polymicrostates
            OUTEEG = pop_getmicrostates(EEG,'Kfrom',nmicrostates,'Kto',nmicrostates,'clustering_algorithm',clustering_algorithms(c),'chronos',0,'MULTI',60,'verbose',0);
            OUTEEGs{c} = OUTEEG;
%            plotsegmentation(OUTEEG);title(num2str(clustering_algorithm));
            %labelings(:,i) = OUTEEG.Z;
            %foo = corr(newerps(:,:,i), OUTEEG.W(:,labelings(:,i))); % corr between frame and every microstate
            %map2frame_correlation(:,i) = diag(foo(:,labelings(:,i))); % mean corr between frame and corresponding microstate
            %map2frame_correlation(:,i) = foo(sub2ind(size(foo),(1:numel(labelings(:,i))).',labelings(:,i)));
            OUTEEG= plotreconstruction(OUTEEG);
            plotactivations(OUTEEG);
            plotmicrostates(OUTEEG);
            recon_error(c,i) = OUTEEG.reconstruction_error;
            %strcat(algoname(c),'_reconstruction')
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            print('-dpng')
%             if OUTEEG.A == 0
%                 recon_error(c,i) = norm(EEG.data - OUTEEG.W*(arr2mat(OUTEEG.Z,OUTEEG.K)),'fro')/norm(EEG.data,'fro');
%             else
%                 recon_error(c,i) = norm(EEG.data - OUTEEG.W*(arr2mat(OUTEEG.Z,OUTEEG.K).*OUTEEG.A),'fro')/norm(EEG.data,'fro');
%             end
        end
    end
end
%%
params.n_iterations = 20;
params.n_init_samples = 20;
params.crit_name = 'cEI';
params.surr_name = 'sStudentTProcessNIG';
params.noise = 1e-6;
params.kernel_name = 'kMaternARD5';
params.kernel_hp_mean = [1];
params.kernel_hp_std = [10];
params.verbose_level = 1;
params.log_filename = 'matbopt.log';
draw=0;
verbose=0;
MULTI = 20;
% beta = 100;
% snr = 10*log10(std(yrec(:))./beta.^(-1/2))
% Y = yrec + rand(J,T)*beta^(-1/2);
max_nits = 1500;
Y = spm_erps(:,:,1);
Y=bsxfun(@rdivide,bsxfun(@minus,Y,mean(Y,2)),std(Y,[],2));
s = 0;
learn_rate_init = 0.2;
learn_rate_decay = 0.001;
save('polypassedparams.mat','Y','K','draw','MULTI','max_nits','verbose','s','learn_rate_init','learn_rate_decay');
% lrateinit, lratedecay, sparsity, smoothness
%F = free_energy_polymicro([0.04,0.0001,-20,.5])
F = free_energy_polymicro([-20,.5])
% 0.5000    0.0001  -10.0000    1.5462
% F = -2.9598
% 0.139327,0.0075485,-7.2494,2.39768
% F = -2.75225
[hyperparams, F] = bayesoptcont('free_energy_polymicro',2,params,[-15 0.01],[0.01 3])
%% Tune variational
params.n_iterations = 20;
params.n_init_samples = 20;
params.crit_name = 'cEI';
params.surr_name = 'sStudentTProcessNIG';
params.noise = 1e-6;
params.kernel_name = 'kMaternARD5';
params.kernel_hp_mean = [1];
params.kernel_hp_std = [10];
params.verbose_level = 1;
params.log_filename = 'matbopt.log';
draw=0;
verbose=0;
MULTI = 5;
beta=1; %hyperparams = 3.7634    0.4341    0.0015
                     % 36.8021    0.5000    0.0003
                     % random data: 64.2085    0.9299    0.4032    0.0001

%beta=100; %hyperparams = 1.0000    0.4127    0.0010 
                        % 468.7576    0.5000    0.0025 
                        % 0.0010    0.2642    0.0029 
                        % 4.86465,0.119712,0.00359363
                        % 3.0214    0.9074    0.4771    0.0010
                        % random data: 29.6228    0.9258    0.1216    0.0012
                        % non-random data : 457.1788    0.9358    0.0441    0.0016
                        %                   393.6398    0.9822    0.4252    0.0041

%beta=1000; %hyperparams = 680.7629    0.5000    0.0028 
                         % 0.4651    0.3923    0.0006
                         % non-random data : 533.6465    0.9000    0.0794    0.0039

%beta = 5000; %hyperparams =  0.0954    0.5000    0.0050
            %               0.5728    0.9002    0.2295    0.0010
%snr = 10*log10(std(yrec(:))./beta.^(-1/2))
%Y = yrec + rand(J,T)*beta^(-1/2);
max_nits = 2000;
kfolds = 1;
alpha = 5;
[~,Z_true] = max(s,[],1);
entropy(Z_true)
K = size(Y,1); %J
save('passedparams.mat','Y','K','draw','alpha','MULTI','G','max_nits','kfolds','verbose','learn_rate_init','learn_rate_decay','Z_true');
%free_energy_variational_wrapper([-3.21099,0.000490851])
[hyperparams, free_energy_test] = bayesoptcont('free_energy_variational_wrapper',2,params,[-10 10^(-10)],[0 0.1])
%% Train both models on grad average
% Zs = zeros(T,nerps);
% Ws = zeros(J,K,nerps);
% As = zeros(K,T,nerps);
nerps = 1;
%OUTEEGs_mono = cell(nerps);
OUTEEGs_poly = cell(nerps);
for i=1:nerps
    EEG.data = mean(spm_erps,3);
    % 2 = kmeans, 3 = agglo, 4 = n-microstates, 5 = dpmeans, 6 = multi-gar
    OUTEEGs_mono{i} = pop_getmicrostates(EEG,'subset',1,'kfrom',K,'kto',K,'clustering_algorithm',6,'multi',40);
    OUTEEGs_poly{i} = pop_getmicrostates(EEG,'subset',1,'kfrom',K,'kto',K,'clustering_algorithm',7,'multi',40);
end
% entropy(OUTEEGs{1}.Z)
% entropy(OUTEEGs{2}.Z)
% mutualinfo(OUTEEGs{1}.Z,OUTEEGs{2}.Z)
%% Mono
[OUTEEGs_mono{1}, com] = merge_correlated_microstates(OUTEEGs_mono{1},4);
plotsegmentation(OUTEEGs_mono{1});
plotactivations(OUTEEGs_mono{1});
plotmicrostates(OUTEEGs_mono{1});
plotmicrostatescorrelation(OUTEEGs_mono{1})
%% Poly
[OUTEEGs_poly{1}, com] = merge_correlated_microstates(OUTEEGs_poly{1});
plotsegmentation(OUTEEGs_poly{1});
plotactivations(OUTEEGs_poly{1});
plotmicrostates(OUTEEGs_poly{1});
plotmicrostatescorrelation(OUTEEGs_poly{1})
