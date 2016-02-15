load('/home/frans/documents/Microstates/EEG_faces/SPMgainmatrix_aceMdspm8_faces_run1_1_c.mat')
% Easy
%J = 128;
% K = 4;

% Hard
J = 32;
K = 7; 
T = 2000;
draw = 1;
%[w,x,s,yrec,p0] = simulate_data(J,K,T,G,draw,EEG.chanlocs);
[w,x,s,yrec,p0] = simulate_data_randomly(J,K,T,draw,EEG.chanlocs); % tilfældige scalp maps, tilfældig tidskomponenter
%%
[~,Z] = max(s,[],1);
minimum_description_length = description_length(Z);
%% 1 run
for i=1:1
MULTI = 10;
verbose = 0;
draw = 0;
alpha = .5;
kfolds = 1;
max_nits = 3000;
beta = 1000;
EEG.data = yrec + rand(J,T)*beta^(-1/2);
gamma2 = 6;
%gamma2 = 0;
learn_rate_init = 0.001;
learn_rate_decay = 0.000001;
K = J;
%[W_var,Mu,M,allZs,beta_var,free_energy,recon_error,m_winner,nits] = variational_microstates_smooth(Y,K,draw,alpha,gamma2,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
OUTEEG_simulated = pop_getmicrostates(EEG,1,K,K,6);
[OUTEEG_simulated, com] = merge_correlated_microstates(OUTEEG_simulated);
plotsegmentation(OUTEEG_simulated);
plotactivations(OUTEEG_simulated);
plotmicrostates(OUTEEG_simulated);
plotmicrostatescorrelation(OUTEEG_simulated)

%plot_true_vs_estimate(w,x,s,W,Mu,M>.5,EEG.chanlocs)
end
%% Compare smoothness to no smoothness
MULTI = 10;
verbose = 0;
draw = 0;
alpha = .5;
kfolds = 1;
max_nits = 3000;
gamma2s = linspace(0,10,2);
betas = logspace(-2,2,2);
nexperiments = 10;
MIs = zeros(size(betas,2),size(gamma2s,2),nexperiments);
[~,Z_true] = max(s,[],1);
for i=1:size(betas,2)
    beta = betas(i); %bad snr
    for j=1:size(gamma2s,2)
        for k=1:nexperiments
            gamma2 = gamma2s(j);
            learn_rate_init = 0.0000001/(exp(abs(gamma2)))
            %learn_rate_init = 0.5;
            learn_rate_decay = 0.00001;
            %gamma2 = log(1 + (p0 - (1-p0)/(K-1))/((1-p0)/(K-1)))-log(1)
            %gamma2 = -5;
            %gamma2 = 0;
            Y = yrec + rand(J,T)*beta^(-1/2);
            %K = size(Y,1);
            K = 7;
            [W_var,Mu,M,allZs,beta_var,free_energy,recon_error,m_winner,nits] = variational_microstates_smooth(Y,K,draw,alpha,gamma2,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
            [~,Z_est] = max(M,[],1);
            Hs(i,j,k) = entropy(Z_est);
            %Ms = arr2mat(Z_est,K);
            MIs(i,j,k) = mutualinfo(Z_true,Z_est);
        end
    end
end
figure;imagesc(mean(MIs,3));colorbar()
xlabel('SNR')
ylabel('Smoothing')
% figure;
% subplot(2,1,1)
% imagesc(s)
% title('True s')
% subplot(2,1,2)
% imagesc(Ms)
%plot_true_vs_estimate(w,x,s,W_var,Mu,M>.5,EEG.chanlocs)
plot(reshape(MIs,5,[])','.')
ylabel('Mutual Info with True Segmentation')
xlabel('Smoothing (gamma_2 = -K) ... No Smoothing (gamma_2 = 0)')
snr = 10*log10(std(yrec(:))./betas.^(-1/2))
meanmis = squeeze(mean(MIs,3));
stdmis = std(MIs,[],3);
figure;boundedline(gamma2s,meanmis,permute(cat(3,meanmis-stdmis,meanmis+stdmis),[1 3 2]),'-')
ylim([2.2 2.6])
legend(['SNR =  ',num2str(snr(1))],['SNR =  ',num2str(snr(2))], ['SNR =  ',num2str(snr(3))], ['SNR =  ',num2str(snr(4))], ['SNR =  ',num2str(snr(5))])
%% Tune hyperparams with BayesOpt, alpha & learning rate & learning rate decay
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
MULTI = 25;
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
snr = 10*log10(std(yrec(:))./beta.^(-1/2))
Y = yrec + rand(J,T)*beta^(-1/2);
max_nits = 1000;
kfolds = 1;
alpha = 0.5;
[~,Z_true] = max(s,[],1);
entropy(Z_true)
learn_rate_decay = 0.0001;
save('passedparams.mat','Y','K','draw','alpha','MULTI','G','max_nits','kfolds','verbose','learn_rate_decay','Z_true');%'learn_rate_init',
%free_energy_variational_wrapper([0.5])
[hyperparams, free_energy_test] = bayesoptcont('free_energy_variational_wrapper',2,params,[0 0.0001],[20 0.1])
%% Performance = f(SNR) %%% LEARNING RATE IS KEY! %%% HYPEROPT FIRST IS VITAL %%%
draw = 0;
verbose = 0;
nb = 5;
lambda = 7;
MULTI = 50;
%MULTI = 2;
nbeta = 5;
%nbeta = 2;
kfolds = 1;
%beta = logspace(-3.5,2.5,nbeta); % whole range, non-random data
beta = logspace(.5,2.5,nbeta); % whole range, random data
%beta = logspace(1.5,2.5,nbeta); % high snr
snr = 10*log10(std(yrec(:))./beta.^(-1/2))
max_nits = 1000;
alpha = 50;
%alpha = 500;
%learn_rate_init = 0.03;
%learn_rate_init = 0.0794;
learn_rate_init = 0.05;
%learn_rate_init = 0.4;
%learn_rate_decay = 0.003;
%learn_rate_decay = 0.0039;
learn_rate_decay = 0.001;
[~,optZ] = max(s);

corrW2w = zeros(nbeta,2);
corrMs2s = zeros(nbeta,2);
corrX2x = zeros(nbeta,2);
recon_err = zeros(nbeta,2);
for b=1:size(beta,2)
    b
    % data
    Y = yrec + rand(J,T)*beta(b)^(-1/2);
    recon_err_denom = norm(yrec,'fro');
    [W_var,Mu,M,beta_var,free_energy,recon_error,m_winner,nits] = variational_microstates_smooth(Y,K,draw,alpha,p0,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
    [~,Ms]=max(M,[],1);
    Z1 = Ms;
    Ms=arr2mat(Ms,K);
    [est_perm1,signs1]=soft_align(W_var,w,Ms,s)
    if any(1-abs(signs1))
        keyboard
    end
    Ms = Ms(est_perm1,:);
    Mu = bsxfun(@times,Mu(est_perm1,:),signs1');
    W_var = bsxfun(@times,W_var(:,est_perm1),signs1);
    corrMs2s(b,1) = mean(abs(diag(corr(Ms',s'))));
    corrW2w(b,1) = mean(abs(diag(corr(W_var,w))));
    for k=1:K
        corrX2x(b,1) = corrX2x(b,1) + abs(corr(x(k,Z1==k)',Mu(k,Z1==k)'));
    end
    %recon_err(b,1) = norm(Y-W_var*(Ms.*Mu),'fro')/recon_err_denom;
    recon_err(b,1) = norm(yrec-W_var*(Ms.*Mu),'fro')/recon_err_denom;
    % pascual-marqui
    %achievable_energy = sum(diag(Y'*Y)-diag((w(:,optZ)'*Y).^2))/(T*(J-1))
    [W,A,Z,~,~] = basicNmicrostates(Y,K,K,MULTI,nb,lambda);
    M2 = arr2mat(Z,K);
    Mu2 = M2.*repmat(A',K,1);
    [est_perm2,signs2]=soft_align(W,w,M2,s)
    if any(1-abs(signs2))
        keyboard
    end
    Mu2 = bsxfun(@times,Mu2(est_perm2,:),signs2');
    M2 = M2(est_perm2,:);
    W = bsxfun(@times,W(:,est_perm2),signs2);
    corrMs2s(b,2) = mean(abs(diag(corr(M2',s'))));
    corrW2w(b,2) = mean(abs(diag(corr(W,w))));
    for k=1:K
        corrX2x(b,2) = corrX2x(b,2) + abs(corr(x(k,Z==k)',A(Z==k)));
    end
    %recon_err(b,2) = norm(Y-W*(Mu2),'fro')/recon_err_denom;
    recon_err(b,2) = norm(yrec-W*(Mu2),'fro')/recon_err_denom;
end
% corrX2x = corrX2x./K;
% figure;
% plot(snr,recon_err(:,1),'k--',snr,corrMs2s(:,1),'r--',snr,corrW2w(:,1),'m--',snr,corrX2x(:,1),'b--')
% ylabel('Mean correlation over K')
% xlabel('SNR in dB')
% legend('Variational Reconstruction Error','Variational Segmentation','Variational Microstates','Variational Activations')
% title('Correlation between estimates and true parameters as a function of SNR')

corrX2x = corrX2x./K;
figure;
plot(snr,recon_err(:,1),'r-',snr,corrMs2s(:,1),'rx',snr,corrW2w(:,1),'ro',snr,corrX2x(:,1),'r.',snr,recon_err(:,2),'b-',snr,corrMs2s(:,2),'bx',snr,corrW2w(:,2),'bo',snr,corrX2x(:,2),'b.')
ylabel('Correlation')
xlabel('SNR in dB')
legend('V Reconstruction Error','V Segmentation','V Microstates','V Activations','N Reconstruction Error','N Segmentation','N Microstates','N Activations')
title('Correlation between estimates and true parameters as a function of SNR')

plot_true_vs_estimate(w,x,s,W_var,Mu,Ms,EEG.chanlocs)
plot_true_vs_estimate(w,x,s,W,Mu2,M2,EEG.chanlocs)

plot_free_energy_and_recon_error(free_energy,recon_error)
%% Gridsearch for lrate
Y = yrec + rand(J,T)*200^(-1/2);
draw = 0;
MULTI = 10;
alpha = 1;
p0 = 0.999;
max_nits = 1000;
kfolds = 4;
% This is how learn_rate is used during training:
% learn_rate = learn_rate_init/(1 + learn_decay);
% M = (1-learn_rate)*oldM+ learn_rate*M;

%learn_rate_inits = 0.1;
learn_rate_inits = logspace(-2,-1,5);
learn_rate_decays = logspace(-4,-3,5);%0.01;%0.00001;
verbose = 0;
best_test_free_energy = zeros(5,5);
for i=1:5
    for j=1:5
        [W_var,Mu,M,beta_var,free_energy,recon_error,m_winner,nits] = variational_microstates_smooth(Y,K,draw,alpha,p0,MULTI,G,max_nits,kfolds,learn_rate_inits(i),learn_rate_decays(j),verbose);
        %[free_energy_train free_energy_test]
        best_test_free_energy(i,j) = mean(free_energy(m_winner,:,2),2);
    end
end

[~,I] = min(min(best_test_free_energy));
bestidx = ind2sub(I,[3,3]);
learn_rate_init = learn_rate_inits(bestidx(1))
learn_rate_decay = learn_rate_decays(bestidx(2))

figure
surf(repmat(learn_rate_inits,5,1),repmat(learn_rate_decays,5,1)',real(best_test_free_energy))
xlabel('Initial lrate')
ylabel('lrate Decay')
zlabel('Test free-energy')
colorbar()

%%

plot_true_vs_estimate(w,x,s,W_var,Mu,M,EEG.chanlocs)

Yrec_hat = W_var*(Mu.*Ms);
figure;
subplot(3,1,1)
plot(yrec',Yrec_hat','.')
xlabel('yrec')
ylabel('yrec_hat')
subplot(3,1,2)
plot(w,W_var(:,est_perm),'.')
xlabel('w')
ylabel('W')
subplot(3,1,3)
plot((x.*s)',(Mu(est_perm,:).*Ms(est_perm,:))','.')
xlabel('x.*s')
ylabel('Mu.*Ms')

figure;
subplot(2,1,1)
plot(1:T,yrec,'r')
title(['Relative reconstruction error is ', num2str(norm(yrec-W_var*(Ms.*Mu),'fro')/norm(yrec,'fro'))])
subplot(2,1,2)
plot(1:T,Yrec_hat,'b')

% subplot(3,1,3)
% plot(1:MULTI,free_energy_pascual)
% y1=get(gca,'ylim');
% hold on
% plot([m_winner m_winner],y1)
% hold off
%% Tune hyperparams with random search, alpha & learning rate
draw=0;
MULTI = 20;
beta=1;
snr = 10*log10(std(yrec(:))./beta.^(-1/2))
Y = yrec + rand(J,T)*beta^(-1/2);
max_nits = 700;
kfolds = 1;
save('passedparams.mat','Y','K','draw','MULTI','G','max_nits','kfolds');
% alpha, p0
%hyperparams = [100 1-eps];
N = 200;
hyperparams = zeros(2,N);
free_energy_test = zeros(1,N);
for i=1:N
    hyperparams(:,i) = rand(2,1).*[0.5,0.5]'+[0.5,0.5]';
    hyperparams(:,i)
    free_energy_test(i) = free_energy_variational_wrapper(hyperparams);
    free_energy_test(i)
end
%epsilon = 10^(-3);
%max_prob = 1-epsilon;
%[hyperparams, free_energy_test] = bayesoptcont('free_energy_variational_wrapper',2,params,[1,2],[0.5,0.9])
scatter(hyperparams(1,:),hyperparams(2,:),1000*ones(N,1),free_energy_test,'.')
colorbar()
%%
b = 5;
lambda = 7;
tic;
[W,A,Z,K,clustering_algorithm] = basicNmicrostates(Y,K,K,MULTI,b,lambda);
toc;
M2 = arr2mat(Z);
X = M2.*repmat(A',K,1);
[est_perm,signs]=soft_align(double(M2>.5),s)
plot_true_vs_estimate(w,x,s,W(:,est_perm),X(est_perm,:),M2(est_perm,:),EEG.chanlocs,signs)
%% Microsoftlost
Y = yrec + rand(J,T)*100^(-1/2);
draw = 1;
Nits = 30;
[W,sig2,p,Nk,loglik] = microsoftlost(Y,K,Nits,draw)

[~,Ms]=max(M,[],1);
Ms=arr2mat(Ms,K);
[est_perm,signs]=soft_align(W_var,w,Ms,s)
plot_true_vs_estimate(w,x,s,W_var(:,est_perm),Mu(est_perm,:),Ms(est_perm,:),EEG.chanlocs(1:(128/J):128),signs)

Yrec_hat = W_var*(Mu.*Ms);
subplot(3,1,1)
plot(yrec',Yrec_hat','.')
xlabel('yrec')
ylabel('yrec_hat')
subplot(3,1,2)
plot(w,W_var(:,est_perm),'.')
xlabel('w')
ylabel('W')
subplot(3,1,3)
plot((x.*s)',(Mu(est_perm,:).*Ms(est_perm,:))','.')
xlabel('x.*s')
ylabel('Mu.*Ms')

figure;
subplot(2,1,1)
plot(1:T,yrec,'r')
title(['Relative reconstruction error is ', num2str(norm(yrec-W_var*(Ms.*Mu),'fro')/norm(yrec,'fro'))])
subplot(2,1,2)
plot(1:T,Yrec_hat,'b')

free_energy_train = mean(free_energy(:,:,1),2);
free_energy_test = mean(free_energy(:,:,2),2);
free_energy_pascual = mean(free_energy(:,:,3),2);

[free_energy_train free_energy_test]

figure;
subplot(2,2,1)
plot(1:MULTI,free_energy_train)
title('Free energy training')
y1=get(gca,'ylim');
hold on
plot([m_winner m_winner],y1)
hold off
subplot(2,2,2)
plot(1:MULTI,free_energy_test)
title('Free energy test')
y1=get(gca,'ylim');
hold on
plot([m_winner m_winner],y1)
hold off
subplot(2,2,3)
plot(1:MULTI,mean(recon_error(:,:,1),2))
title('Relative reconstruction error training')
y1=get(gca,'ylim');
hold on
plot([m_winner m_winner],y1)
hold off
subplot(2,2,4)
plot(1:MULTI,mean(recon_error(:,:,2),2))
title('Relative reconstruction error test')
y1=get(gca,'ylim');
hold on
plot([m_winner m_winner],y1)
hold off
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