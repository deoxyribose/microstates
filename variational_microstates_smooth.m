function  [W,Mu,M,sigma2,allZs,beta,free_energy,recon_error,m_winner,nits,varargin] = variational_microstates_smooth(Yall,varargin)
switch nargin
    case 11
        K = varargin{1};
        draw = varargin{2};
        alpha = varargin{3};
        gamma2 = varargin{4};
        MULTI = varargin{5};
        max_nits = varargin{6};
        trigger = varargin{7};
        learn_rate_init = varargin{8};
        learn_rate_decay = varargin{9};
        verbose = varargin{10};
    otherwise
end
[allJ,T]=size(Yall);
normYall = norm(Yall,'fro');
converged = 0;
nits = 0;
best_free_energy = Inf;
free_energy = zeros(MULTI,1);
recon_error = zeros(MULTI,1);
% train multi models and evaluate each on a test set, select the model that
% has lowest mean test-error over folds
bestM = rand(K,T);
bestW = zeros(allJ,K);
count=1;
Nits = 10000;
if draw==1,
    h = figure;
    t_arr=zeros(1,ceil(Nits/10));
    beta_arr=zeros(1,ceil(Nits/10));
    lrate_arr=zeros(1,ceil(Nits/10));
    mact=zeros(K,ceil(Nits/10));
    xact=zeros(K,ceil(Nits/10));
    err=zeros(1,ceil(Nits/10));
    no_Mu=zeros(K,ceil(Nits/10));
    no_W=zeros(K,ceil(Nits/10));
end
nits_total = 0;
disp('Fitting model...')
minimum_description_length = Inf;
% test & training
kfolds = 1;
splitj = floor(allJ/kfolds);
jidxs = randperm(allJ);
jidxs = reshape(jidxs,splitj,kfolds);
allZs = zeros(T,MULTI);
for m=1:MULTI
    %learn_rate0 = max(learn_rate_init + (rand-.5)*0.5,0.1);
    learn_rate0 = learn_rate_init;
    learn_decay0 = learn_rate_decay;
    %     learn_rate0 = learn_rate_init;
    %     learn_decay0 = learn_rate_decay;
    if verbose
        disp(['Running Variational Microstates with ', num2str(K),' clusters for the ', num2str(m), 'th time'])
    end
    random_stuff_1 = rand(K,T)>0.5;
    %
    %
    %     %M = Ms+random_stuff_1(:,1:t);
    %     M = random_stuff_1+rand(K,T)*0.3;
    %     M = bsxfun(@rdivide,M,sum(M,1));
    for k=1:1
        learn_rate = learn_rate0;
        learn_decay = learn_decay0;
        nits = 0;
        converged = 0;
        if kfolds == 1
            train_jidx = 1:allJ;
            test_jidx = [];
        else
            train_jidx = jidxs(:,[1:kfolds]~=k);
            test_jidx = jidxs(:,k);
        end
        Y = Yall(train_jidx,:);
        [J,T]=size(Y);
        %%% INIT %%%
        [~,~,U]=svd(Y',0);
        W=U(:,1:K);
        W=bsxfun(@times,W,1./sqrt(sum(W.^2,1)));
        Mu = sqrt((Y'*W).^2);
        Mu = Mu';
        [~,Ms] = max(Mu,[],1);
        Ms = arr2mat(Ms,K);
        M = Ms+rand(K,T)*0.1;
        M = bsxfun(@rdivide,M,sum(M,1));

        beta=1/mean(mean(Y.*Y));

        WcircW = W.*W;
        diagWW=sum(WcircW,1);
        WY = W'*Y;
        sigma2 = 1./((1/alpha^2)+beta*bsxfun(@times,Ms,diagWW'));
        normY = norm(Y,'fro');
        recon_err_train = norm(Y-W*(Ms.*Mu),'fro')/normY;
        while ~converged
            % Latent/Z/S/M
            oldM = M;
            oldMs = Ms;
            %             M=softmax(beta*( (WY).*Mu - 0.5*diag(diagWW)*(Mu.^2 + sigma2))+ gamma1 + gamma2*([M(:,2:end) zeros(K,1)] + [zeros(K,1) M(:,1:end-1)]));
            Mtm1 = [zeros(K,1) M(:,1:end-1)];
            Mtp1 = [M(:,2:end) zeros(K,1)];
            %gamma2 = (1-learn_rate)*gamma2 -learn_rate*sum(sum(Mtm1.*M));
            %gamma2 = 1/mean(mean(Mtm1.*M));
            M=softmax(beta*( (WY).*Mu - 0.5*diag(diagWW)*(Mu.^2 + sigma2))+gamma2*(Mtm1+Mtp1));
            if isnan(M)
                keyboard
            end
            
            % find lowest
            
            % We expect a new microstate when the change in Mu
            % is large, but we'll add noise too.
            %             derivMoX = diff(Mu.*Ms,1,2);
            %             D = abs(derivMoX);
            %             D(:,end+1) = D(:,end);
            %             D = D + mean(mean(abs(derivMoX)));
            %             D = bsxfun(@times,D,1./(max(D,[],2)+0.1));
            %D = rand(K,T)>0.1;
            D = 1;
            % How seriously we take the new developement depends on the
            % learning rate, which is large in the beginning, small towards
            % convergence
            M = (1-learn_rate*D).*oldM + learn_rate*D.*M;
            
            [~,Z]=max(M,[],1);
            oldMs = Ms;
            Ms=arr2mat(Z,K);
            
            % Microstates/W/Gamma
            oldW = W;
            
            %             D2 = log(cumsum(Ms,2));
            %             D2(D2<0) = 0.1;
            %             D2 = repmat(mean(D2,2)',J,1);
            %             D2 = D2 + (rand(J,K)<0.1)*0.1;
            %             D2 = bsxfun(@times,D2,1./(max(D2,[],2)+0.1));
            %D2 = rand(J,K)>0.2;
            D2 = 1;
            W = bsxfun(@rdivide,Y*(Mu.*M)',sum((Mu.^2+sigma2).*M,2)');
            if isnan(W)
                keyboard
            end
            %W = (1-learn_rate)*D2.*oldW + learn_rate*D2.*W;
            %
            %W=bsxfun(@rdivide,bsxfun(@minus,W,mean(W,1)),std(W,[],1)); % normalize
            WcircW = W.*W;
            diagWW=sum(WcircW,1);
            WY = W'*Y;
            denom = (1/alpha^2)+beta*bsxfun(@times,M,diagWW');
            %sigma2 = learn_rate*sigma2+ learn_rate*1./denom;
            sigma2 = 1./denom;
            if isnan(sigma2)
                keyboard
            end
            % Activations/Mu/X/A
            oldMu = Mu;
            Mu = (1-learn_rate)*Mu+ learn_rate*((beta*(WY.*M))./denom);
            if isnan(Mu)
                keyboard
            end
            % beta
            oldbeta = beta;
            beta = J*T/sum(sum(Y.^2-2*Y.*(W*(Mu.*M))+WcircW*((Mu.^2 + sigma2).*M)));
            if isnan(beta)
                keyboard
            end
            Yrec=W*(Mu.*Ms);
            recon_err_train_old = recon_err_train;
            recon_err_train = norm(Y-Yrec,'fro')/normY;
            
            if (draw==1)&&(rem(nits,10)==1),
                t_arr(count)=nits_total;
                no_W(count)=norm(W,'fro');
                no_Mu(count)=norm(Mu,'fro');
                beta_arr(count)=beta;
                lrate_arr(count) = learn_rate;
                err(count)=norm(Y-Yrec,'fro')/norm(Y,'fro');
                mact(:,count)=sum(Ms,2);
                xact(:,count)=sum(Mu,2);
                
                figure(h)
                subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
                subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
                %subplot(3,2,3),plot(t_arr(1:count),no_Mu(1:count)), title('||Mu||')
                subplot(3,2,3),plot(t_arr(1:count),xact(:,1:count)'),title('||Mu_k||')
                %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
                grid
                subplot(3,2,4),plot(t_arr(1:count),mact(:,1:count)'),title('||M_k||')
                subplot(3,2,5), plot(t_arr(1:count),lrate_arr(1:count)),title('lrate')
                subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count),'k'),title('\beta')
                hold on
                drawnow
                count=count+1;
            end
            learn_rate = learn_rate0/(1 + learn_decay)^nits;
            %%% CONVERGENCE CHECK
            %converged_allK = all(sum(Ms,2)>0);
            %             converged_beta = abs(beta-oldbeta)/oldbeta<tol;
            %             noldM = norm(oldM);
            %             converged_M = abs(norm(M)-noldM)/noldM<tol;
            %             converged = converged_M*converged_allK*converged_allT;
            %             converged = converged_M*converged_allK*converged_beta;
            %converged =  abs(recon_err_train_old-recon_err_train) < tol;
            tol = 10^(-3);
            converged = abs(beta-oldbeta)/oldbeta< tol;
%             if nits>2
%                 if isnan(free_energy(nits,m)) | nits > max_nits
%                     if verbose
%                         disp('Number of iterations exceeds limit, proceeding without convergence')
%                     end
%                     converged = 1;
%                 else
%                     converged = norm(M-oldM)/norm(M) < tol;
%                 end
%             end
            %converged = norm(M-oldM)/norm(M) < tol;
            if converged*verbose
                disp(['Converged in ', num2str(nits), ' iterations']);
            end
            if nits > max_nits
                if verbose
                    disp('Number of iterations exceeds limit, proceeding without convergence')
                end
                converged = 1;
            end
            nits = nits +1;
            nits_total = nits_total +1;
        end
        
        % TRAINING FREE ENERGY
        sum_K_T = sum(sum(-1/2*log(2*pi*sigma2) + M.*log(M) + (Mu.^2 + sigma2)/(2*alpha^2)-gamma2*M.*Mtm1));
        sum_J_T = (beta/2)*sum(sum(Y.^2 - 2*Y.*(W*(Mu.*M)) + W.^2*((Mu.^2 + sigma2).*M)));
        constants = T/2*(K*(log(2*pi*alpha^2*K)-1)-J*log(beta/(2*pi))); %- gamma1*KT;
        free_energy(m) = sum_K_T + sum_J_T + constants;
        recon_error(m) = norm(Y-W*(Ms.*Mu),'fro')/norm(Y,'fro');
        
        %         % TEST FREE ENERGY
        %         Y_test = Yall(test_jidx,:);
        %         %W_pred = Sigma12invSigma22*W;
        %         W_pred = bsxfun(@rdivide,Y_test*(Mu.*M)',sum((Mu.^2+sigma2).*M,2)');
        %         WcircW_pred = W_pred.*W_pred;
        %         diagWW_pred=sum(WcircW_pred,1);
        %         WY_pred = W_pred'*Y_test;
        %         beta_pred = beta;
        %         M_pred = M;
        %         Mu_pred = Mu;
        %         sigma2_pred = sigma2;
        %         J_test = numel(test_jidx);
        %         W_pred = bsxfun(@rdivide,Y_test*(Mu_pred.*M_pred)',sum((Mu_pred.^2+sigma2_pred).*M_pred,2)');
        
        %         sum_K_T_test = sum(sum(-1/2*log(2*pi*sigma2_pred) + M_pred.*log(M_pred) + (Mu_pred.^2 + sigma2_pred)/(2*alpha^2)));
        %         sum_J_T_test = (beta_pred/2)*sum(sum(Y_test.^2 - 2*Y_test.*(W_pred*(Mu_pred.*M_pred)) + W_pred.^2*((Mu_pred.^2 + sigma2_pred).*M_pred)));
        %         constants_test = T/2*(K*(log(2*pi*alpha^2*K)-1)-J_test*log(beta_pred/(2*pi)));
        %         free_energy(m,k,2) = sum_K_T_test + sum_J_T_test + constants_test;
        %         recon_error(m,k,2) = norm(Y_test-W_pred*(Ms.*Mu),'fro')/norm(Y_test,'fro');
        
        % PASCUAL FREE ENERGY
        %         Wtmp=bsxfun(@times,W,1./sqrt(sum(W.^2,1)));
        %         free_energy(m,k,3) = sum(diag(Y'*Y)-diag((Wtmp(:,Z)'*Y).^2))/(T*(J-1));
        %         if draw
        %             y1=get(gca,'ylim');
        %             hold on
        %             plot([nits_total nits_total],y1,'k')
        %             drawnow
        %         end
    end
    if free_energy(m) < best_free_energy
    %encoding_length = description_length(Z);
    %if encoding_length < minimum_description_length;
        %recon_error(m) < best_free_energy
        %minimum_description_length = encoding_length;
        best_free_energy = free_energy(m);
        bestMu = Mu;
        bestM = M;
        bestsigma2 = sigma2;
        bestW = bsxfun(@rdivide,Yall*(bestMu.*bestM)',sum((bestMu.^2+bestsigma2).*bestM,2)');
        m_winner = m;
    else %if best_free_energy == Inf
          %  keyboard
       % end
    end
    allZs(:,m) = Z;
end
W = bestW;
Mu = bestMu;
M = bestM;
beta = allJ*T/sum(sum(Yall.^2-2*Yall.*(W*(Mu.*M))+(W.*W)*((Mu.^2 + bestsigma2).*M)));
% disp(num2str(mean(free_energy(m_winner,:,3))));'        sum_K_T = sum(sum(-1/2*log(2*pi*sigma2) + M.*log(M) + (Mu.^2 + sigma2)/(2*alpha^2)));
sum_J_T = (beta/2)*sum(sum(Yall.^2 - 2*Yall.*(W*(Mu.*M)) + W.^2*((Mu.^2 + sigma2).*M)));
constants = T/2*(K*(log(2*pi*alpha^2*K)-1)-J*log(beta/(2*pi)));
sum_K_T + sum_J_T + constants
disp('Done')
end