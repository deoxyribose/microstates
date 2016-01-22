function  [W,X,M,beta,gamma1,gamma2,nits] = polymicro_smooth(Yall,K,draw,gamma1,gamma2)
[Jall,T]=size(Yall);
splitj = floor(Jall/4);
MULTI = 1;
best_test_err = 1;
X=rand(K,T);
Y = Yall;
for i=1:MULTI
     i
%     jidxs = randperm(Jall);
%     test_jidx = jidxs(1:splitj);
%     train_jidx = jidxs(splitj+1:end);
%     Y_test = Yall(test_jidx,:);
%     Y = Yall(train_jidx,:);
%     coords_train = chanlocs(train_jidx);
%     coords_test = chanlocs(test_jidx);
%     xcoord_train = [coords_train.X]';
%     ycoord_train = [coords_train.Y]';
%     zcoord_train = [coords_train.Z]';
%     xcoord_test = [coords_test.X]';
%     ycoord_test = [coords_test.Y]';
%     zcoord_test = [coords_test.Z]';
    
    [J,T]=size(Y);
    [~,~,U]=svd(Y',0);
    W=U(:,1:K);
    M=0.1*rand(K,T);
    M(M<0.01)=0.01;
    beta=1/mean(mean(Y.*Y));
    count=1;
    Nits = 3000;
    test_err=ones(1,ceil(Nits/10));
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
    nits = 1;
    
    logistic = @(x)1/(exp(-x)+1);
    while ~converged
        
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
        for t = 1:size(X,2)
            tmp = Wsq*diag(M(:,t)) + diag((1-M(:,t)).*diagww);%+0.00001*eye(K);
            if rcond(tmp) < 10e-10
                keyboard
            end
%             D = abs(randi(2,K,1)-1);
            %D = abs(rand(K,1));
%             X(:,t) = D.*X(:,t) + (1-D).*(tmp\(WY(:,t)));
            X(:,t) = (tmp\(WY(:,t)));
        end
        
        %%% W
        %     tmp = zeros(K,1);
        %     for k = 1:K
        %         for t = 1:T
        %             tmp(k) = tmp(k) + X(k,t)^2*(M(k,t) - M(k,t)^2);
        %         end
        %     end
        MoX = M.*X;
        
        % %    W = Y*MoX'/((MoX*MoX') + MoX*((1-M).*X)');
        %     for t = 1:size(X,2)
        %         W(:,t) = Y(:,t)*MoX/((MoX*MoX') + diag(MoX*((1-M).*x)')));
        %     end
        if rcond(MoX*MoX') < 10e-5
            disp('Ill-conditioned MoX*MoX!')
        end
        
        %W = 0.95*W + 0.05*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
        
        %%% dropout
%         D = abs(randi(2,J,K)-1);
        %     W = D.*W + (1-D).*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
        
        %%% soft dropout
        %D = abs(rand(J,K));
        %     for i=1:K
        %         most_correlated_maps = abs(corr(W(:,1),Y))>0.5;
        %         W(:,i) = Y(:,most_correlated_maps)*MoX(i,most_correlated_maps)'/(MoX(i,most_correlated_maps)*MoX(i,most_correlated_maps)' + diag(sum(X.^2.*(M - M.^2),2)));
        %         W(:,i) = D.*oldW(:,i) + (1-D).*W(:,i);
        %     end
%         W = D.*W + (1-D).*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
 W = (Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
        
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
        
%         D = abs(randi(2,K,T)-1);
        %D = repmat(D,1,T);
        %M = D.*M + (1-D).*arrayfun(logistic,beta*X.*(bsxfun(@times,X.*(M-.5),diagww) + W'*(Y-yrec)) + log((gamma10*gamma01)/gamma00^2) + ([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)])*log((gamma00*gamma11)/(gamma10*gamma01)));
        M = 0.95*M + 0.05*arrayfun(logistic,beta*X.*(-0.5*bsxfun(@times,X.*(1-2*M),diag(W'*W)) + W'*(Y-yrec)) + gamma1 + gamma2*([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)]));
        
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
%         D = abs(rand(J,T));
%         yrec = D.*(W*(M.*X)) + (1-D).*yrec;
        yrec = W*(M.*X);
        beta_old=beta;
        beta = J*T/(sum(sum((Y-yrec).^2)) + sum(sum(diag(diag(W'*W))*(X.^2).*(M-M.^2))));
        
        
        
        if ( abs(beta-beta_old)/beta_old )< 10^(-5)*beta_old,
            nits
            break
        end
%         if rem(nits,10)==1
%             W_test = zeros(size(test_jidx,2),K);
%             for k=1:K
%                 F = scatteredInterpolant(xcoord_train,ycoord_train,zcoord_train,W(:,k));
%                 W_test(:,k) = F(xcoord_test,ycoord_test,zcoord_test);
%             end
%             yrec_test = W_test*(M.*X);
%             test_err(count)=norm(Y_test-yrec_test)/norm(Y_test);
%             %             if count>1
%             %                 converged = abs(test_err(count)-test_err(count)) < 0.001;
%             %             end
%             count=count+1;
%         end
        if (draw==1)&&(rem(nits,10)==1),
            Ms=M>0.5;
            t_arr(count)=nits;
            no_W(count)=norm(W);
            no_X(count)=norm(X);
            beta_arr(count)=beta;
            err(count)=norm(Y-yrec)/norm(Y);
            mact(:,count)=sum(Ms,2);
            xact(:,count)=sum(X,2);
            subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
            subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
            %subplot(3,2,3),plot(t_arr(1:count),no_X(1:count)), title('||Mu||')
            subplot(3,2,3),plot(xact(:,1:count)')
            %subplot(3,2,4),plot(t_arr(1:count),test_err(1:count)),title('Test Rel Rec Error')
            %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
            grid
            subplot(3,2,5),plot(mact(:,1:count)')
            subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count)),title('\beta')
            drawnow
            count=count+1;
        end
        %     tol = 10^(-5);
        %     converged1 = norm(M - oldM)/norm(oldM) < tol;
        %     converged2 = norm(W - oldW)/norm(oldW) < tol;
        %     converged3 = norm(X - oldX)/norm(oldX) < tol;
        %     converged = converged1*converged2*converged3;
        nits = nits +1;
        if nits>1500
            converged = 1;
        end
    end
end
%     end
%     if test_err(count) < best_test_err
%         best_X = X;
%         best_M = M;
%         best_W = zeros(size(Yall,1),K);
%         for k=1:K
%             F = scatteredInterpolant(xcoord_train,ycoord_train,zcoord_train,W(:,k));
%             best_W(:,k) = F([chanlocs.X]',[chanlocs.Y]',[chanlocs.Z]');
%         end
%         best_test_err = test_err(count)
%     end
% end
% converged = 0;
% count=1;
% draw = 1;
% nits=1;
% Y = Yall;
% J = Jall;
% X = best_X;
% M = best_M;
% W = best_W;
% 
% if draw==1,
%     t_arr=zeros(1,ceil(Nits/10));
%     beta_arr=zeros(1,ceil(Nits/10));
%     mact=zeros(K,ceil(Nits/10));
%     err=zeros(1,ceil(Nits/10));
%     no_X=zeros(K,ceil(Nits/10));
%     no_W=zeros(K,ceil(Nits/10));
%     figure(1)
% end
% 
% while ~converged
%     Wsq = W'*W;
%     diagww = diag(Wsq);
%     WY = W'*Y;
%     D = randi(2,K,1)-1;
%     for t = 1:size(X,2)
%         tmp = (Wsq*diag(M(:,t)) + diag((1-M(:,t)).*diagww));
%         if rcond(tmp) < 10e-10
%             keyboard
%         end
%         X(:,t) = D.*X(:,t) + (1-D).*(tmp\(WY(:,t)));
%     end
%     MoX = M.*X;
%     D = abs(randi(2,J,K)-1);
%     W = D.*W + (1-D).*(Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
%     Wsq = W'*W;
%     diagww = diag(Wsq);
%     WY = W'*Y;
%     oldX = X;
%     D = randi(2,K,1)-1;
%     for t = 1:size(X,2)
%         tmp = (Wsq*diag(M(:,t)) + diag((1-M(:,t)).*diagww));
%         X(:,t) = D.*X(:,t) + (1-D).*(tmp\(WY(:,t)));
%     end
%     yrec = W*(M.*X);
%     oldM = M;
%     D = abs(randi(2,K,T)-1);
%     M = D.*M + (1-D).*arrayfun(logistic,beta*X.*(-0.5*bsxfun(@times,X.*(1-2*M),diag(W'*W)) + W'*(Y-yrec)) + gamma1 + gamma2*([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)]));
%     D = abs(rand(J,T));
%     yrec = D.*(W*(M.*X)) + (1-D).*yrec;
%     
%     beta_old=beta;
%     beta = J*T/(sum(sum((Y-yrec).^2)) + sum(sum(diag(diag(W'*W))*(X.^2).*(M-M.^2))));
%     %     if (draw==1)&&(rem(nits,10)==1),
%     %         Ms=M>0.5;
%     %         t_arr(count)=nits;
%     %         no_W(count)=norm(W);
%     %         no_X(count)=norm(X);
%     %         beta_arr(count)=beta;
%     %         err(count)=norm(Y-yrec)/norm(Y);
%     %         mact(:,count)=sum(Ms,2);
%     %         xact(:,count)=sum(X,2);
%     %         subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
%     %         subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
%     %         %subplot(3,2,3),plot(t_arr(1:count),no_X(1:count)), title('||Mu||')
%     %         subplot(3,2,3),plot(xact(:,1:count)')
%     %         %subplot(3,2,4),plot(t_arr(1:count),test_err(1:count)),title('Test Rel Rec Error')
%     %         %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
%     %         grid
%     %         subplot(3,2,5),plot(mact(:,1:count)')
%     %         subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count)),title('\beta')
%     %         drawnow
%     %         if count>2
%     %             converged = abs(err(count)-err(count)) < 0.1;
%     %         end
%     %         count=count+1;
%     %     end
%     
%     nits = nits +1;
%     
%     if ( abs(beta-beta_old)/beta_old )< 10^(-5)*beta_old, break, end
% end
% title('DONE!');
% end