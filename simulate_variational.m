J = 128;
K = 7;
T = 1000;
draw = 0;
[w,x,s,yrec,p0] = simulate_data(J,K,T,draw);
% 
% % Variational should find this minimum:
% sigma2 = eps*ones(K,T); alpha = eps;
% s(s==0) = eps;
% sum_K_T = sum(sum(-1/2*log(2*pi*sigma2) + s.*log(s) + (x.^2 + sigma2)/(2*alpha^2)));
% sum_J_T = beta/2*sum(sum(Y.^2 - 2*Y.*(w*(x.*s)) + w.^2*((x.^2 + sigma2).*s)));
% constants = T/2*(K*(log(2*pi*alpha^2*K)-1)-J*log(beta/(2*pi)));
% free_energy = sum_K_T + sum_J_T + constants
% 
% % Basic should find this minimum:
% [~,Z] = max(s,[],1);
% sum(diag(Y'*Y)-diag((w(:,Z)'*Y).^2))/(T*(J-1))
%% Performance = f(SNR)
draw = 0;
nb = 5;
lambda = 7;
MULTI = 50;
nbeta = 50;
beta = logspace(-1,3,nbeta);
snr = 10*log10(std(yrec(:))./beta.^(-1/2));
%corrX2x = zeros(10,2);
corrW2w = zeros(nbeta,2);
corrMs2s = zeros(nbeta,2);
corrX2x = zeros(nbeta,2);
for b=1:size(beta,2)
    Y = yrec + rand(J,T)*beta(b)^(-1/2);
    [W_var,Mu,M,~,~,~] = variational_microstates_smooth(Y,K,draw,0.8,0.999,MULTI);
    [~,Ms]=max(M,[],1);
    Ms=arr2mat(Ms,K);
    [est_perm,signs]=soft_align(Ms,s)
    Ms = Ms(est_perm,:);
    Mu = Mu(est_perm,:);
    W_var = W_var(:,est_perm);
    %plot_true_vs_estimate(w,x,s,W_var,Mu,Ms,[1:J],[1:T],EEG.chanlocs)
    corrMs2s(b,1) = mean(abs(diag(corr(Ms',s'))));
    corrW2w(b,1) = mean(abs(diag(corr(W_var,w))));
    for k=1:K
        corrX2x(b,1) = corrX2x(b,1) + abs(corr(x(k,Z1==k)',Mu(k,Z1==k)'));
    end
    [W,A,Z,~,~] = basicNmicrostates(Y,K,K,MULTI,nb,lambda);
    M2 = arr2mat(Z);
    mewtwo = M2.*repmat(A',K,1);
    mewtwo = mewtwo(est_perm,:);
    [est_perm,signs]=soft_align(M2,s)
    M2 = M2(est_perm,:);
    W = W(:,est_perm);
    corrMs2s(b,2) = mean(diag(corr(M2',s')));
    corrW2w(b,2) = mean(abs(diag(corr(W,w))));
    for k=1:K
        corrX2x(b,2) = corrX2x(b,2) + abs(corr(x(k,Z==k)',A(Z==k)));
    end
end
corrX2x = corrX2x./K; 
figure;
plot(snr,corrMs2s(:,1),snr,corrMs2s(:,2),snr,corrW2w(:,1),snr,corrW2w(:,2),snr,corrX2x(:,1),snr,corrX2x(:,2))
ylabel('Mean correlation over K')
xlabel('SNR in dB')
legend('Variational Segmentation','N-microstates Segmentation','Variational Microstates','N-microstates')
title('Correlation between estimates and true parameters as a function of SNR')
%% Compare basic N-microstates with Variational Microstates
Y = yrec + rand(J,T)*100^(-1/2);
draw = 0;
b = 5;
lambda = 7;
MULTI = 50;
tic;
[W_var,Mu,M,beta,free_energy,nits] = variational_microstates_smooth(Y,K,draw,0.8,0.999,MULTI);
toc;
[~,Ms]=max(M,[],1);
Ms=arr2mat(Ms,K);
[est_perm,signs]=soft_align(Ms,s)
plot_true_vs_estimate(w,x,s,W_var(:,est_perm),Mu(est_perm,:),Ms(est_perm,:),[1:J],[1:T],EEG.chanlocs)

tic;
[W,A,Z,K,clustering_algorithm] = basicNmicrostates(Y,K,K,MULTI,b,lambda);
toc;
M2 = arr2mat(Z);
X = M2.*repmat(A',K,1);
[est_perm,signs]=soft_align(double(M2>.5),s)
plot_true_vs_estimate(w,x,s,W(:,est_perm),X(est_perm,:),M2(est_perm,:),[1:J],[1:T],EEG.chanlocs)
%% Training and test
% CV on electrodes because we need the chronology intact for smoothing
splitj = floor(J/4);
jidxs = randperm(J);
% splitt = floor(T/10);
% tidxs = randperm(T);
test_jidx = jidxs(1:splitj);
% test_tidx = tidxs(1:splitt);
train_jidx = jidxs(splitj+1:end);
% train_tidx = tidxs(splitt+1:end);
w_test = w(test_jidx,:);
% Y_test = Y(test_jidx,test_tidx);
% Y_train = Y(train_jidx,train_tidx);
Y_test = Y(test_jidx,:);
Y_train = Y(train_jidx,:);


Sigma12invSigma22 = (G(test_jidx,:)*G(train_jidx,:)')*inv(G(train_jidx,:)*G(train_jidx,:)');
%% best achievable mse
W_pred = Sigma12invSigma22*w(train_jidx,:);
mse_true = mean(mean((W_pred-w_test).^2))
%% initializing model randomly
draw = 1;
[W,Mu,M,beta,free_energy,nits] = variational_microstates_smooth(Y_train,K,draw,0.8,p0);
W_pred = Sigma12invSigma22*W;
mse = mean(mean((W_pred-w_test).^2))
plot_true_vs_estimate(w,x,s,W,X,M,train_jidx,[1:T],EEG.chanlocs)
%% BayesOpt
params.n_iterations = 10;
params.n_init_samples = 10;
params.crit_name = 'cEI';
params.surr_name = 'sStudentTProcessNIG';
params.noise = 1e-6;
params.kernel_name = 'kMaternARD5';
params.kernel_hp_mean = [1];
params.kernel_hp_std = [10];
params.verbose_level = 1;
params.log_filename = 'matbopt.log';
draw=0;
save('passedparams.mat','Y_train','K','draw','Sigma12invSigma22','w_test');
%mse = free_energy_polymicro_wrapper2(gammas);
[gammas, mse] = bayesoptcont('free_energy_polymicro_wrapper2',2,params,[-10,0],[0,4])

[W,X,M,Beta,gamma1,gamma2,nits] = polymicro_smooth(Y_train,K,draw,gammas(1),gammas(2),rand(K,T));
plot_true_vs_estimate(w,x,s,W,X,M,train_jidx,[1:T],EEG.chanlocs)
%% Example
EEG = pop_loadset('filename','simulated3.set','filepath','/home/frans/MATLAB packages/eeglab13_4_4b/plugins/microstates0.1/');
%load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat')
%EEG.chanlocs = chanlocs(1:30);
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5); % auto-clean plugin
EEG=pop_chanedit(EEG, 'load',{'/home/frans/documents/Microstates/bdf-2-M-Neuch/biosemi128.xyz' 'filetype' 'autodetect'});
pop_chanedit(gcbf, [], 'return', []);
EEG.data = Y;
EEG.pnts = size(Y,2);
eeglab('redraw')

[~,true_idx] = max(s,[],1);

min_set_size = 10;
models = [6];
[training_size, error_train, error_test, map_generalization, map_potential] = CVv2(EEG, Y, models, min_set_size, true_idx, w);
figure; 
plot(training_size,mean(error_train,3),'c',training_size,mean(error_test,3),'m',training_size,mean(map_generalization,3),'y', training_size,mean(map_potential,3),'k')
xlabel('Training set size')
ylabel('Error')
legend('Training error','Test error','Proportion of mislabeled microstates using estimated maps','Proportion of mislabeled microstates using true maps')
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
