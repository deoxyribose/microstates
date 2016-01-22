function  [W,X,M,beta,nits] = polymicro(Y,K,draw)
%[W,Mu,M,beta] = multi_garrote(Y,K,Nits,alpha,draW)
%
[J,T]=size(Y);
[~,~,U]=svd(Y',0);
W=U(:,1:K);
M=(1/K)*ones(K,T); %M=0.45 + 0.1*rand(K,T);

%%% Pass true M
% s = zeros(K,T);
% s(1,:) = sin((1:T)/10)>0.8;
% %s(2,:) = sin((1:T)/100)>0.5;
% s(3,:) = sin((1:T)/1000)>0.2;
% s(2,find(sum(s,1) == 0)) = 1;
% M = s;
% oldM = M;

X = rand(K,T);
%%% Pass true X
% X = 0.1*bsxfun(@times,sin(bsxfun(@plus,bsxfun(@times,repmat(1:T,K,1),[1,2,3]'/10^2),[5,10,15]')) ...
%     + cos(bsxfun(@plus,bsxfun(@times,repmat(1:T,K,1),[3,2,1]'/10^2),[2,4,19]')),[17,4,5]');
% oldX = X;

beta=1/mean(mean(Y.*Y));
gamma = -0.52;

count=1;
Nits = 3000;
if draw==1,
    t_arr=zeros(1,ceil(Nits/10));
    beta_arr=zeros(1,ceil(Nits/10));
    mact=zeros(K,ceil(Nits/10));
    err=zeros(1,ceil(Nits/10));
    no_X=zeros(K,ceil(Nits/10));
    no_W=zeros(K,ceil(Nits/10));
    figure(1)
end
converged = 0;
nits = 0;

logistic = @(x)1/(exp(-x)+1);
while ~converged
    
    %%% W
%     tmp = zeros(K,1);
%     for k = 1:K
%         for t = 1:T
%             tmp(k) = tmp(k) + X(k,t)^2*(M(k,t) - M(k,t)^2);
%         end
%     end
    MoX = M.*X;
    oldW = W;
    
% %    W = Y*MoX'/((MoX*MoX') + MoX*((1-M).*X)');
%     for t = 1:size(X,2)
%         W(:,t) = Y(:,t)*MoX/((MoX*MoX') + diag(MoX*((1-M).*x)')));
%     end
    if rcond(MoX*MoX') < 10e-10
        keyboard
    end
    
    %W = 0.95*W + 0.05*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
    
    %%% dropout
%     D = abs(randi(2,J,K)-1);
%     W = D.*W + (1-D).*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
    
    %%% soft dropout
     D = abs(rand(J,K));
%     for i=1:K
%         most_correlated_maps = abs(corr(W(:,1),Y))>0.5;
%         W(:,i) = Y(:,most_correlated_maps)*MoX(i,most_correlated_maps)'/(MoX(i,most_correlated_maps)*MoX(i,most_correlated_maps)' + diag(sum(X.^2.*(M - M.^2),2)));
%         W(:,i) = D.*oldW(:,i) + (1-D).*W(:,i);
%     end
    W = D.*W + (1-D).*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
    
%%% X
%     led1 = zeros(K,T);
%     led2 = zeros(K,T);
%     for t = 1:T
%         wsq = zeros(K,1);
%         for l = 1:K
%             for j = 1:J
%                 for k = 1:K
%                     led1(l,t) = led1(l,t) + W(j,k)*W(j,l)*M(k,t)*X(k,t);
%                     end
%                 wsq(l) = wsq(l) + W(j,l)^2;
%             end
%             led2(l,t) = X(l,t)*(1-M(l,t))*wsq(l);
%         end
%     end
% 
    Wsq = W'*W;
    diagww = diag(Wsq);
    WY = W'*Y;
    oldX = X;
    for t = 1:size(X,2)
        tmp = (Wsq*diag(M(:,t)) + diag((1-M(:,t)).*diagww));
        if rcond(tmp) < 10e-10
            keyboard
        end
        %D = abs(randi(2,K,1)-1);
        D = abs(rand(K,1));
        X(:,t) = D.*X(:,t) + (1-D).*(tmp\(WY(:,t)));
    end
    %%% M
    yrec = W*(M.*X);
    oldM = M;
%     for l=1:K
%         for t=1:T
%             led1 = 0;
%             led2 = 0;
%             for j=1:J
%                 led1 = led1 + (Y(j,t)-yrec(j,t))*W(j,l)*X(l,t);
%                 led2 = led2 + W(j,l)^2*X(l,t)^2*(1-2*M(l,t));
%             end
%             M(l,t) = 0.95*M(l,t) + 0.05*logistic(gamma+beta*led1-beta/2*led2);
%         end
%     end
    
    %M = 0.95*M + 0.05*arrayfun(logistic,gamma + beta*X.*(bsxfun(@times,X.*(M-.5),diagww) + W'*(Y-yrec)));
    %D = abs(randi(2,K,T)-1);
    D = abs(rand(K,T));
    M = D.*M + (1-D).*arrayfun(logistic,gamma + beta*X.*(bsxfun(@times,X.*(M-.5),diagww) + W'*(Y-yrec)));
    
    %%% beta
%     
%     led1 = 0;
%     led2 = 0;
%     for t = 1:T        
%         for j = 1:J
%             yrec = 0;
%             for k = 1:K
%                 yrec = yrec + M(k,t)*W(j,k)*X(k,t);
%                 led1 = led1 + W(j,k)^2*X(k,t)^2*M(k,t)^2-M(k,t);
%             end
%             led2 = led2 + (Y(j,t)-yrec);
%         end
%     end

    %yrec = W*(M.*X);

    %D = abs(randi(2,J,T)-1);
    D = abs(rand(J,T));
    yrec = D.*(W*(M.*X)) + (1-D).*yrec;
    beta = 0.95*beta + 0.05*J*T/(sum(sum((Y-yrec).^2)) - sum(sum(W.^2*(X.^2.*(M.^2-M)))));
    
    sumsumM = sum(sum(M));
    gamma = 0.95*gamma - 0.05*log(sumsumM/(K*T-sumsumM));
    
    if (draw==1)&&(rem(nits,10)==1),
        Ms=M>0.5;
        t_arr(count)=nits;
        no_W(count)=norm(W);
        no_X(count)=norm(X);
        beta_arr(count)=beta;
        gamma_arr(count)=gamma;
        err(count)=norm(Y-yrec)/norm(Y);
        mact(:,count)=sum(Ms,2);
        
        subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
        subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
        subplot(3,2,3),plot(t_arr(1:count),no_X(1:count)), title('||Mu||')
        %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
        subplot(3,2,5),plot(t_arr(1:count),gamma_arr(1:count)), title('\gamma')
        grid
        subplot(3,2,4),plot(mact(:,1:count)')
        subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count)),title('\beta')
        drawnow
        count=count+1;
    end
    tol = 10^(-4);
    converged1 = norm(M - oldM)/norm(oldM) < tol;
    converged2 = norm(W - oldW)/norm(oldW) < tol;
    converged3 = norm(X - oldX)/norm(oldX) < tol;
    converged = converged1*converged2*converged3;
    nits = nits +1;
    if nits>1500
        converged = 1;
    end
end