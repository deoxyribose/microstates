function  [W,Mu,M,beta,nits] = multi_garrote(Y,K,alpha,draw,uni)
%[W,Mu,M,beta] = multi_garrote(Y,K,Nits,alpha,draw)
%
[d,N]=size(Y);
[~,~,U]=svd(Y',0);
W=U(:,1:K);
M=(1/K)*ones(K,N);  %0.45+0.1*rand(K,N);
beta=1/mean(mean(Y.*Y));

%
count=1;
Nits = 1000;
if draw==1,
    t_arr=zeros(1,ceil(Nits/10));
    beta_arr=zeros(1,ceil(Nits/10));
    mact=zeros(K,ceil(Nits/10));
    err=zeros(1,ceil(Nits/10));
    no_Mu=zeros(K,ceil(Nits/10));
    no_W=zeros(K,ceil(Nits/10));
end
[J, K] = size(W);
T = size(Y,2);
Mu = zeros(K,T);
converged = 0;
nits = 0;
while ~converged
    %disp(['t = ',int2str(t),' of ',int2str(Nits)])
    WcircW = W.*W;
    diagWW=sum(WcircW,1);
    WY = W'*Y;
    %     Mu2 = zeros(K,T);
    %     for k=1:K
    %         for t=1:T
    %             numerator_sum = 0;
    %             denominator_sum = 0;
    %             for j=1:J
    %                 numerator_sum = numerator_sum + Y(j,t)*W(j,k)*M(k,t);
    %                 denominator_sum = denominator_sum + (W(j,k)^2)*M(k,t);
    %             end
    %             Mu2(k,t) = (beta*numerator_sum)/(1/alpha^2 + beta*denominator_sum);
    %         end
    %     end
    sigma2 = 1./((1/alpha^2)+beta*bsxfun(@times,M,diagWW'));
    oldMu = Mu;
    Mu = beta*(WY.*M).*sigma2;
    
    %     denominator = 0;
    %     for j = 1:J
    %         for t = 1:T
    %             first_sum = 0;
    %             second_sum = 0;
    %             for k = 1:K
    %                 first_sum = first_sum + W(j,k)*Mu(k,t)*M(k,t);
    %                 second_sum = second_sum + (W(j,k)^2)*(Mu(k,t)^2 + sigma2(k,t))*M(k,t);
    %             end
    %             denominator = denominator + Y(j,t)^2 - 2*Y(j,t)*first_sum+second_sum;
    %         end
    %     end
    %     beta2 = (J*T)/denominator;
    beta=0.9*beta+0.1*(d*N)/sum(sum(Y.*Y- 2*Y.*(W*(Mu.*M)) +WcircW*((Mu.^2 + sigma2).*M)));% a kind of line-search?
    %     M2 = zeros(K,T);
    %     for t = 1:T
    %         for k = 1:K
    %             notanothersum = 0;
    %             for j = 1:J
    %                 notanothersum = notanothersum + Y(j,t)*W(j,k)*Mu(k,t)-0.5*(W(j,k)^2)*(Mu(k,t)^2 + sigma2(k,t));
    %             end
    %             M2(k,t) = exp(-beta*notanothersum);
    %         end
    %         M2(:,t) = M2(:,t)./sum(M2(:,t));
    %     end
    oldM = M;
    M=0.95*M+ 0.05*softmax(uni*beta*( (WY).*Mu - 0.5*diag(diagWW)*(Mu.^2 + sigma2)));
    %
    oldW = W;
    W = bsxfun(@rdivide,Y*(Mu.*M)',sum((Mu.^2+sigma2).*M,2)');
    %     W2 = zeros(J,K);
    %     for j=1:J
    %         for k=1:K
    %             numerator = 0;
    %             denominator = 0;
    %             for t=1:T
    %                 numerator = numerator + Y(j,t)*Mu(k,t)*M(k,t);
    %                 denominator = denominator + (Mu(k,t)^2 + sigma2(k,t))*M(k,t);
    %             end
    %             W2(j,k) =  numerator/denominator;
    %         end
    %     end
    [~,Ms]=max(M,[],1);
    Ms=arr2mat(Ms,K);
    Yrec=W*(Mu.*Ms);
    if (draw==1)&&(rem(nits,10)==1),
        t_arr(count)=nits;
        no_W(count)=norm(W);
        no_Mu(count)=norm(Mu);
        beta_arr(count)=beta;
        err(count)=norm(Y-Yrec)/norm(Y);
        mact(:,count)=sum(Ms,2);
        
        %figure(1)
        subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
        subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
        subplot(3,2,3),plot(t_arr(1:count),no_Mu(1:count)), title('||Mu||')
        %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
        grid
        subplot(3,2,4),plot(mact(:,1:count)')
        subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count)),title('\beta')
        drawnow
        count=count+1;
    end
    tol = 10^(-3);
    converged1 = all(all((M - oldM) < tol))
    converged2 = all(all((W - oldW) < tol))
    converged3 = all(all((Mu - oldMu) < tol))
    converged = converged1*converged2*converged3;
    %converged = (std(sum(Ms,2)) < 26) & (nits > 50);
%    converged = nits > 15;
    nits = nits +1
%     if nits==100
%         keyboard
%         figure
%         imagesc(Ms)
%     end
end













