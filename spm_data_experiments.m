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
n = 10;
[J,T,nepochs] = size(Y);
% divide Y into cells of JxTxn, and one cell with the number of epochs
% remaining after division by n, and mean over epochs
spm_erps = cell2mat(cellfun(@(x)mean(x,3),mat2cell(Y, [J],[T],[ones(1,floor(nepochs/n))*n rem(nepochs,n)]),'UniformOutput',false));
clustering_algorithm = 7;
nerps = size(spm_erps,3);
%% Tune hyperparams for polymicrostates
Y = mean(spm_erps,3);
Y=bsxfun(@rdivide,bsxfun(@minus,Y,mean(Y,2)),std(Y,[],2));

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
    %OUTEEGs_mono{i} = pop_getmicrostates(EEG,1,K,K,6);
    OUTEEGs_poly{i} = pop_getmicrostates(EEG,1,K,K,7);
end
% entropy(OUTEEGs{1}.Z)
% entropy(OUTEEGs{2}.Z)
% mutualinfo(OUTEEGs{1}.Z,OUTEEGs{2}.Z)
%% Mono
[OUTEEGs_mono{1}, com] = merge_correlated_microstates(OUTEEGs_mono{1});
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
