 load('/home/frans/documents/Microstates/EEG_faces/SPMgainmatrix_aceMdspm8_faces_run1_1_c.mat')
% Easy
% J = 128;
% K = 4;

% Hard
J = 32;
K = 7; 
T = 1000;
draw = 1;
%[w,x,s,yrec,p0] = simulate_data(J,K,T,G,draw,EEG.chanlocs);
[w,x,s,yrec,p0] = simulate_data_randomly(J,K,T,draw,EEG.chanlocs); % tilfældige scalp maps, tilfældig tidskomponenter
[~,Z_true] = max(s,[],1);
%% Convergence plot
beta = 1000;
MULTI = 10;
verbose = 1;
draw = 1;
alpha = .5;
kfolds = 1;
max_nits = 1500;
learn_rate_init = 0.2;
%learn_rate_init = 0.5;
learn_rate_decay = 0.001;
Y = yrec + rand(J,T)*beta^(-1/2);
[W_var,Mu,M,beta_var,free_energy,recon_error,m_winner,nits] = variational_microstates_smooth(Y,K,draw,alpha,p0,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
[~,Z_est] = max(M,[],1);
Ms = arr2mat(Z_est,K);
plot_true_vs_estimate(w,x,s,W_var,Mu,Ms,EEG.chanlocs)
%% Model selection
beta = 1000;
MULTI = 10;
verbose = 1;
draw = 1;
alpha = .5;
kfolds = 1;
max_nits = 1500;
learn_rate_init = 0.2;
%learn_rate_init = 0.5;
learn_rate_decay = 0.001;
Y = yrec + rand(J,T)*beta^(-1/2);
for k=4:10
    [W_var,Mu,M,beta_var,free_energy,recon_error,m_winner,nits] = variational_microstates_smooth(Y,K,draw,alpha,p0,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
end

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
MULTI = 30;
%beta=1; %hyperparams = 3.7634    0.4341    0.0015
                     % 36.8021    0.5000    0.0003
                     % random data: 64.2085    0.9299    0.4032    0.0001

beta=100; %hyperparams = 1.0000    0.4127    0.0010 
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
save('passedparams.mat','Y','K','draw','MULTI','G','max_nits','kfolds','verbose','Z_true');
%free_energy_variational_wrapper([0.5 0.01 0.0001])
[hyperparams, free_energy_test] = bayesoptcont('free_energy_variational_wrapper',4,params,[0 0.00001],[0 0.1])
%% Performance = f(SNR)
draw = 0;
verbose = 0;
nb = 3;
nk = 3;
nalgo = 3;
lambda = 6.5;
MULTI = 25;
%MULTI = 2;
nbeta = 50;
%nbeta = 2;
kfolds = 1;
%beta = logspace(-3.5,2.5,nbeta); % whole range, non-random data
beta = logspace(.5,3,nbeta); % whole range, random data
%beta = logspace(1.5,2.5,nbeta); % high snr
snr = 10*log10(std(yrec(:))./beta.^(-1/2))
max_nits = 1000;
alpha = 0.5;
learn_rate_init = 0.01;
learn_rate_decay = 0.0001;
[~,optZ] = max(s);
n = 100; % number of experiments
free_energy = zeros(nbeta,nk,nalgo,MULTI,n);
recon_error = zeros(nbeta,nk,nalgo,MULTI,n);
MIs = zeros(nbeta,nk,nalgo,MULTI,n);
Hs = zeros(nbeta,nk,nalgo,MULTI,n);
winners = zeros(nbeta,nk,nalgo,n);
ks = [4,7,10];
%k = 7;
for b=1:size(beta,2)
    for k=1:3
        for e = 1:n
            Y = yrec + rand(J,T)*beta(b)^(-1/2);
            [~,~,~,allZs_V,beta_V,free_energy(b,k,1,:,e),recon_error(b,k,1,:,e),winners(b,k,1,e),nits] = variational_microstates_smooth(Y,ks(k),draw,alpha,0,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
            MIs(b,k,1,:,e) = cellfun(@(x)mutualinfo(Z_true,x),num2cell(allZs_V,1));
            Hs(b,k,1,:,e) = cellfun(@(x)entropy(x),num2cell(allZs_V,1));
            [~,~,~,allZs_Vs,beta_V,free_energy(b,k,2,:,e),recon_error(b,k,2,:,e),winners(b,k,2,e),nits] = variational_microstates_smooth(Y,ks(k),draw,alpha,lambda,MULTI,G,max_nits,kfolds,learn_rate_init,learn_rate_decay,verbose);
            MIs(b,k,2,:,e) = cellfun(@(x)mutualinfo(Z_true,x),num2cell(allZs_Vs,1));
            Hs(b,k,2,:,e) = cellfun(@(x)entropy(x),num2cell(allZs_V,1));
            [~,~,~,allZs_N,free_energy(b,k,3,:,e),recon_error(b,k,3,:,e),winners(b,k,3,e)] = basicNmicrostates(Y,ks(k),ks(k),MULTI,nb,lambda);
            MIs(b,k,3,:,e) = cellfun(@(x)mutualinfo(Z_true,x),num2cell(allZs_N,1));
            Hs(b,k,3,:,e) = cellfun(@(x)entropy(x),num2cell(allZs_N,1));
        end
    end
end
mutualinfoZZ = mutualinfo(Z_true,Z_true);
%mutualinfo_random_Z = mutualinfo(Z_true,randi(K,size(Z_true)));

% size(MIs) = nbeta x nk x nalgorithms x nruns x nexperiments
% size(winners) = nbeta x nk x nalgorithms
for b=1:nbeta
     for k=1:3
         for a=1:3
             for e=1:n
                 if 
                MIs_of_winners(b,k,a,e)=MIs(b,k,a,winners(b,k,a,e),e);
             end
         end
    end
end
% sub1 = repmat(1:nbeta,3*2);
% sub2 = repmat(1:3,nbeta*2);
% ind = sub2ind( size(MIs), sub1(:), [1 2 1 2], [1 3 2 2]);
% K = resize(J(ind), [2 2]);

% figure;plot(squeeze(free_energy(5,2,1,:)),squeeze(MIs(5,2,1,:)),'*')
% figure;plot(squeeze(free_energy(5,2,2,:)),squeeze(MIs(5,2,2,:)),'*')
% figure;plot(squeeze(recon_error(5,2,1,:)),squeeze(MIs(5,2,1,:)),'*')
% figure;plot(squeeze(recon_error(5,2,2,:)),squeeze(MIs(5,2,2,:)),'*')

mean_experimental_result = mean(MIs_of_winners,4); % nbeta x nk x nalgorithms
var_experimental_result = var(MIs_of_winners,1,4); 
std_error = 3*var_experimental_result./sqrt(n); % nbeta x nk x nalgorithms
errorbars = permute(std_error,[1 4 2 3]);

%errorbars = permute(cat(4,abs(quantile(MIs,.25,4)-MIs_of_winners), abs(quantile(MIs,.75,4)-MIs_of_winners)),[1 4 2 3]); % nbeta x nsides x nlines x nalgorithms
figure;
ax1 = subplot(1,3,1)
boundedline(snr,mean_experimental_result(:,:,1),errorbars(:,:,:,1),'r-')
axis([0 10 0 3])
ylabel('I(Z_{true},Z_{est}) in bits')
xlabel('SNR in dB')
ax2 = subplot(1,3,2)
boundedline(snr,mean_experimental_result(:,:,2),errorbars(:,:,:,2),'g-')
%axis([ax1 ax2],[0 10 0 3])
axis([0 10 0 3])
ax2 = subplot(1,3,3)
boundedline(snr,mean_experimental_result(:,:,3),errorbars(:,:,:,3),'b-')
%axis([ax1 ax2],[0 10 0 3])
axis([0 10 0 3])
%line([0, 10],[mutualinfoZZ,mutualinfoZZ],'Color','g','Linewidth',2)
ylabel('I(Z_{true},Z_{est}) in bits')
xlabel('SNR in dB')
legend('','V K=4','','V K=7','','V K=10','','N K=4','','N K=7','','N K=10','H(Z_{true})')
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