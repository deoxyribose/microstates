function  [W,X,M,allZs,beta,free_energy,recon_error,m_winner,gamma1,gamma2,nits] = polymicro_smooth(Y,varargin)
for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~isstr(Param)
        error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'k'
            K = Value;
        case 'draw'
            draw = Value;
        case 'multi'
            MULTI = Value;
        case 'max_nits'
            max_nits = Value;
        case 'learn_rate_init'
            learn_rate_init = Value;
        case 'learn_rate_decay'
            learn_rate_decay = Value;
        case 'gamma1'
            gamma1 = double(Value);
        case 'gamma2'
            gamma2 = double(Value);
        case 'verbose'
            verbose = Value;
        otherwise
            error(['Unknown input parameter ''' Param ''' ???'])
    end
end
[Jall,T]=size(Y);
splitj = floor(Jall/4);
X=rand(K,T);
free_energy = zeros(max_nits,MULTI);
free_energy_end = zeros(1,MULTI);
recon_error_end = zeros(1,MULTI);
recon_error = zeros(max_nits,MULTI);
allZs = zeros(T,MULTI);
%best_free_energy = -Inf;
minimum_description_length = Inf;
Nits = 3000;
if draw==1,
    h = figure;
    t_arr=zeros(1,ceil(Nits/10));
    beta_arr=zeros(1,ceil(Nits/10));
    lrate=zeros(1,ceil(Nits/10));
    mact=zeros(K,ceil(Nits/10));
    err=zeros(1,ceil(Nits/10));
    no_X=zeros(K,ceil(Nits/10));
    no_W=zeros(K,ceil(Nits/10));
    figure(1)
end
for m=1:MULTI
    learn_rate0 = max(learn_rate_init + (rand-.5)*0.5,0.1);
    learn_rate = learn_rate0;
    learn_decay = learn_rate_decay;
    [J,T]=size(Y);
    [~,~,U]=svd(Y',0);
    W=U(:,1:K);
    M=0.1*rand(K,T);
    M(M<0.01)=0.01;
    beta=1/mean(mean(Y.*Y));
    count=1;
    test_err=ones(1,ceil(Nits/10));
    converged = 0;
    nits = 1;
    
    %logistic = @(x)1/(exp(-x)+1);
    while ~converged
        
        Wsq = W'*W;
        diagww = diag(Wsq);
        WY = W'*Y;
        for t = 1:size(X,2)
            X(:,t) = ((Wsq*diag(M(:,t)) + diag((1-M(:,t)).*diagww))\(WY(:,t)));
        end
        
        MoX = M.*X;
        %         if rcond(MoX*MoX') < 10e-5
        %             disp('Ill-conditioned MoX*MoX!')
        %         end
        W = (Y*MoX'/(MoX*MoX' + diag(sum(X.^2.*(M - M.^2),2))));
        
        %%% M
        yrec = W*(M.*X);
        oldM = M;
        %M = (1-learn_rate)*M + learn_rate*arrayfun(logistic,beta*X.*(-0.5*bsxfun(@times,X.*(1-2*M),diag(W'*W)) + W'*(Y-yrec)) + gamma1 + gamma2*([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)]));
        M = (1-learn_rate)*M + learn_rate./(exp(-(beta*X.*(-0.5*bsxfun(@times,X.*(1-2*M),diag(W'*W)) + W'*(Y-yrec)) + gamma1 + gamma2*([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)])))+1);
        yrec = W*(M.*X);
        beta_old=beta;
        beta = J*T/(sum(sum((Y-yrec).^2)) + sum(sum(diag(diag(W'*W))*(X.^2).*(M-M.^2))));
        
        
        
        learn_rate = learn_rate0/(1 + learn_decay)^nits;
        tol = 10^(-2);
        %converged_M = all(sum(M,2)>0.01)*all(sum(M,1)>0.01);
        converged_beta = ( abs(beta-beta_old)/beta_old )< tol*beta_old;
        if nits>2
            if isnan(free_energy(nits,m)) | nits > max_nits
                if verbose
                    disp('Number of iterations exceeds limit, proceeding without convergence')
                end
                converged = 1;
            else
                %converged_recon_err =  abs(recon_error(nits-1,m)-recon_error(nits,m)) < tol;
                %converged_free_energy =  abs(free_energy(nits-1,m)-free_energy(nits,m))/abs(free_energy(nits-1,m)) < tol;
                %converged = converged_free_energy*converged_beta;
                converged = norm(M-oldM)/norm(M) < tol;
            end
        end
        if converged*verbose
            disp(['Converged in ', num2str(nits), ' iterations']);
        end
        if (draw==1)&&(rem(nits,10)==1),
            Ms=M>0.5;
            t_arr(count)=nits;
            no_W(count)=norm(W);
            no_X(count)=norm(X);
            beta_arr(count)=beta;
            lrate(count)=learn_rate;
            err(count)=norm(Y-yrec)/norm(Y);
            mact(:,count)=sum(Ms,2);
            xact(:,count)=sum(X,2);
            figure(h)
            subplot(3,2,1),plot(t_arr(1:count),no_W(1:count)),title('||W||')
            %subplot(3,2,2),plot(t_arr(1:count),err(1:count)),title('Rel Rec Error')
            subplot(3,2,2),plot(1:nits-1,free_energy(1:nits-1,m)),title('Free energy')
            %subplot(3,2,3),plot(t_arr(1:count),no_X(1:count)), title('||Mu||')
            subplot(3,2,3),plot(xact(:,1:count)'),title('||Mu_k||')
            subplot(3,2,4),plot(t_arr(1:count),lrate(1:count)),title('Learning Rate')
            %subplot(3,2,5),plot3(U(:,1)'*Y,U(:,2)'*Y,U(:,3)'*Y,'.')
            grid
            subplot(3,2,5),plot(mact(:,1:count)'), title('||M_k||')
            subplot(3,2,6), plot(t_arr(1:count),beta_arr(1:count)),title('\beta')
            drawnow
            count=count+1;
        end

        %gamma1 = log((gamma10*gamma01)/gamma00^2)
              % = log(gamma10*gamma01)-log(gamma00^2)
              % = log(gamma10)+log(gamma01)-log(gamma00^2)
              % = log(1-gamma11)+log(1-gamma00)-log(gamma00^2)
        % gamma1 = log((gamma10*gamma01)/gamma00^2)
               % = log((gamma10/gamma00)*(gamma01/gamma00))
               % = log(gamma10/gamma00) + log(gamma01/gamma00)
               % = log((1-gamma11)/gamma00) + log((1-gamma00)/gamma00)
        %exp(gamma1) = (1-gamma11-gamma00+gamma00*gamma11)/gamma00^2
        %gamma2 = log((gamma00*gamma11)/(gamma10*gamma01))
              % = log((gamma00*gamma11)/((1-gamma11)*(1-gamma00)))
              % = log((gamma00*gamma11)/(1-gamma11-gamma00+gamma00*gamma11))
        %exp(gamma2) = (gamma00*gamma11)/(1-gamma11-gamma00+gamma00*gamma11)

        % from sympy import *
        % expgamma1,expgamma2,gamma00,gamma11 = symbols('expgamma1 expgamma2 gamma00 gamma11')
        % solve([Eq((1-gamma11-gamma00+gamma00*gamma11)/gamma00^2,expgamma1),Eq((gamma00*gamma11)/(1-gamma11-gamma00+gamma00*gamma11),expgamma2)],[gamma00, gamma11])

        expgamma1 = exp(gamma1);
        expgamma2 = exp(gamma2);
        gamma00 = (expgamma1*expgamma2^2 + expgamma2 - sqrt(expgamma2^2*(expgamma1^2*expgamma2^2 - 2*expgamma1*expgamma2 + 4*expgamma1 + 1)))/(2*expgamma1*expgamma2*(expgamma2 - 1));
        gamma11 = (expgamma1*expgamma2^2 + expgamma2 - sqrt(expgamma2^2*(expgamma1^2*expgamma2^2 - 2*expgamma1*expgamma2 + 4*expgamma1 + 1)))/(2*(expgamma2 - 1));
        gamma01 = max(1-gamma00,0);
        gamma10 = min(1-gamma11,1);
        F1 = .5*J*T*log(beta/(2*pi))-beta/2*sum(sum((Y-W*(X.*M)).^2))-beta/2*sum(sum(diag(diag(W'*W))*(X.^2.*(M-M.^2))));
        F2 = log(gamma00) + sum(sum(M*log(gamma10/gamma00) + [M(:,end) M(:,1:end-1)]*log(gamma01/gamma00) + M.*[M(:,end) M(:,1:end-1)]*log((gamma00*gamma11)/(gamma01*gamma10))));
        F3 = - sum(sum(M.*log(M) + (1-M).*log(1-M)));
        
        free_energy(nits,m) = double(F1+F2+F3);
        free_energy_end(m) = free_energy(nits,m);
        %allZs(:,m) = Z;
        nits = nits+1;
        recon_error(nits,m) = norm(Y-yrec,'fro')/norm(Y,'fro');
        recon_error_end(m) = recon_error(nits,m);
    end
    Z = M>.5;
    encoding_length = description_length(Z);
    if encoding_length < minimum_description_length;
        %recon_error(m) < best_free_energy
        %free_energy_end(m) > best_free_energy
        minimum_description_length = encoding_length;
        %best_free_energy = free_energy_end(m);
        bestM = M;
        m_winner = m;
    end
end
end

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

%%% W
%     tmp = zeros(K,1);
%     for k = 1:K
%         for t = 1:T
%             tmp(k) = tmp(k) + X(k,t)^2*(M(k,t) - M(k,t)^2);
%         end
%     end


% %    W = Y*MoX'/((MoX*MoX') + MoX*((1-M).*X)');
%     for t = 1:size(X,2)
%         W(:,t) = Y(:,t)*MoX/((MoX*MoX') + diag(MoX*((1-M).*x)')));
%     end

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