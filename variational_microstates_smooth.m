function  [W,Mu,M,beta,free_energy,nits] = variational_microstates_smooth(Y,K,draw,alpha,p0,MULTI)

[J,T]=size(Y);

%
count=1;
Nits = 3000;
if draw==1,
    t_arr=zeros(1,ceil(Nits/10));
    beta_arr=zeros(1,ceil(Nits/10));
    mact=zeros(K,ceil(Nits/10));
    err=zeros(1,ceil(Nits/10));
    no_Mu=zeros(K,ceil(Nits/10));
    no_W=zeros(K,ceil(Nits/10));
end
converged = 0;
nits = 0;
best_free_energy = Inf;
gamma1 = log((p0*K - 1)/(K-1));
gamma2 = log((1-p0)/(K-1));
for m=1:MULTI
converged = 0;
nits = 0;

% INITIALIZATION
[~,~,U]=svd(Y',0);
W=U(:,1:K);
Mu = sqrt((Y'*W).^2);
Mu = Mu';
% M = zeros(K,T);
% for k=1:K
%     tmp = sort(Mu(:,k));
%     M(k,:) = Mu(:,k) > tmp(5);
% end
% M = bsxfun(@rdivide,M,sum(M));
[~,s] = max(Mu,[],1);
M = arr2mat(s,K)+rand(K,T);
M = bsxfun(@rdivide,M,sum(M));
% best_init_energy_so_far = Inf;
% energy = 0;
% for m=1:3
% [W,lbl,~,energy] = kmeans_fast(Y',K);
%     if energy < best_init_energy_so_far
%         bestWkm = W;
%         bestlbl = lbl;
%         best_init_energy_so_far = energy;
%     end
% end
% W = bestWkm' + rand(J,K)*0.1;
% W = bsxfun(@minus,bsxfun(@rdivide,W,std(W,[],1)),mean(W,1));
% 
% 

%W = datasample(Y,K,2); % take K random frames
%W = bsxfun(@rdivide,W,sqrt(sum(W.^2,1)));

sigma2 = rand(K,T)/10;
beta=1/mean(mean(Y.*Y));



while ~converged

        WcircW = W.*W;
    diagWW=sum(WcircW,1);
    WY = W'*Y;

        % Latent/Z/S/M
    oldM = M;
    M=0.95*M+ 0.05*softmax(beta*( (WY).*Mu - 0.5*diag(diagWW)*(Mu.^2 + sigma2))+ gamma1 + gamma2*([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)]));

    [~,Ms]=max(M,[],1);
    Ms=arr2mat(Ms,K);
    
    % Microstates/W/Gamma
    oldW = W;
    W = bsxfun(@rdivide,Y*(Mu.*M)',sum((Mu.^2+sigma2).*M,2)');
    WcircW = W.*W;
    diagWW=sum(WcircW,1);
    WY = W'*Y;

    sigma2 = 1./((1/alpha^2)+beta*bsxfun(@times,M,diagWW'));
    
    % Activations/Mu/X/A
    oldMu = Mu;
    Mu = beta*(WY.*M).*sigma2;
    
    % beta
    oldbeta = beta;
    beta=J*T/sum(sum(Y.*Y- 2*Y.*(W*(Mu.*M)) +WcircW*((Mu.^2 + sigma2).*M)));
        
    Yrec=W*(Mu.*Ms);
    if (draw==1)&&(rem(nits,10)==1),
        t_arr(count)=nits;
        no_W(count)=norm(W);
        no_Mu(count)=norm(Mu);
        beta_arr(count)=beta;
        err(count)=norm(Y-Yrec)/norm(Y);
        mact(:,count)=sum(Ms,2);
        xact(:,count)=sum(Mu,2);
        
        %figure(1)
        subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
        subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
        %subplot(3,2,3),plot(t_arr(1:count),no_Mu(1:count)), title('||Mu||')
        subplot(3,2,3),plot(xact(:,1:count)')
        %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
        grid
        subplot(3,2,4),plot(mact(:,1:count)')
        subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count)),title('\beta')
        drawnow
        count=count+1;
    end
     tol = 10^(-4);
%     converged1 = all(all((M - oldM) < tol));
%     converged2 = all(all((W - oldW) < tol));
%     converged3 = all(all((Mu - oldMu) < tol));
%     converged = converged1*converged2*converged3;
    converged_allK = all(sum(Ms,2)>0);
    converged_beta = abs(beta-oldbeta)/oldbeta<tol;
    converged = converged_beta*converged_allK;
    %converged = (std(sum(Ms,2)) < 26) & (nits > 50);
%    converged = nits > 15;
     nits = nits +1;
%     if nits==100
%         keyboard
%         figure
%         imagesc(Ms)
%     end
end
sum_K_T = sum(sum(-1/2*log(2*pi*sigma2) + M.*log(M) + (Mu.^2 + sigma2)/(2*alpha^2)));
sum_J_T = beta/2*sum(sum(Y.^2 - 2*Y.*(W*(Mu.*M)) + W.^2*((Mu.^2 + sigma2).*M)));
constants = T/2*(K*(log(2*pi*alpha^2*K)-1)-J*log(beta/(2*pi)));
free_energy = sum_K_T + sum_J_T + constants;
if free_energy < best_free_energy
    best_free_energy = free_energy
    bestW = W;
    bestMu = Mu;
    bestM = M;
end
end
W = bestW;
Mu = bestMu;
M = bestM;
end