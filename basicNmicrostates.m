function [W,A,Z,allZs,smu2,recon_err,winner,K,clustering_algorithm,scv,sgcv,mcv] = basicNmicrostates(Y,kfrom,kto,MULTI,b,lambda,eps,verbose)
clustering_algorithm = 'N-microstates';
[J,T] = size(Y);
mcv = zeros(kto-kfrom+1,1);
scv = zeros(kto-kfrom+1,1);
smu2 = zeros(kto-kfrom+1,MULTI);
best_mcv_so_far = Inf;
S = Y*Y';
[eigvecs,eigvals] = eig(S);
c=1;
for k=kfrom:kto
    sgcv(c) = trace(eigvals(k+1:J-1,k+1:J-1))/(J-1);
    sgcv(c) = sgcv(c)/((J-1-k)/(J-1))^2;
    c=c+1;
end
c=1;
Yall = Y;
for k=kfrom:kto
    best_energy_so_far = Inf;
    allZs = zeros(T,MULTI);
    cv_finished = 1;
    for m=1:MULTI
        if kto ~= kfrom
            cv_finished = 0;
            j = 1;
            Y_j = Yall(j,:);
            Y = Yall([1:j-1 j+1:J],:);
        end
        if verbose
            disp(['Running N-microstates with ', num2str(k),' clusters for the ', num2str(m), 'th time'])
        end
        % The basic N-microstate algorithm
        s2 = 0;
        %eps = 10^(-5);
        % init
        W = datasample(Y,k,2); % take random frames
        W = bsxfun(@rdivide,W,sqrt(sum(W.^2,1)));
        % step 3
        [~,Z] = max((Y'*W).^2,[],2);
        allK = 0;
        smu2(k-kfrom+1,m) = sum(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(T*(J-1));
        %         smu22 = 0;
        %         for t=1:T
        %             smu22 = smu22 + (Y(:,t)'*Y(:,t) - (W(:,Z(t))'*Y(:,t)).^2);
        %         end
        %         smu22 = smu22/(T*(J-1));
        while abs(s2 - smu2(k-kfrom+1,m))>eps*smu2(k-kfrom+1,m) & cv_finished %allK % minimize orthogonal squared distance between each observed vector and corresponding microstates
            s2 = smu2(k-kfrom+1,m);
            [~,Z] = max((Y'*W).^2,[],2);
            W = zeros(J,k);
            %S = zeros(J,J,k);
            for i = 1:k;
                %S(:,:,i) = Y(:,Z==i)*Y(:,Z==i)'; % sum covariances of all spatial maps for microstate i
                [W(:,i),~] = eigs(Y(:,Z==i)*Y(:,Z==i)',1); % In what direction is the variance within the microstate largest? Eigs returns normalized eigenvector corresp. to largest eigenvalue
            end
            % W = bsxfun(@rdivide,W,sqrt(sum(W.^2,1)));
            smu2(k-kfrom+1,m) = sum(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(T*(J-1));
            %scv2 = norm(,2)^2;
            %allK = numel(unique(Z)) == k;
        end
        % find activations as projections of observed vectors onto microstates
        %         a = zeros(k,J);
        %         for t=1:T;
        %             kappa = Z(t);
        %             for i=1:k
        %                 if kappa == i
        %                     a(kappa,t) = Y(:,t)'*W(:,kappa);
        %                 end
        %             end
        %         end
        A = diag(Y'*W(:,Z));
        recon_err(m) = norm(Y-W*(arr2mat(Z',k).*repmat(A',k,1)),'fro')/norm(Y,'fro');
        if kto ~= kfrom
            scv(c) = scv(c) + norm(Y_j-W(j,Z)*A,2)^2/J;
            if j+1<J
                j = j+1;
            else
                cv_finished = 1;
            end
        end
        sD2 = sum(diag(Y'*Y))/(T*(J-1));
        R2 = 1 - smu2(k-kfrom+1,m)/sD2;
        allZs(:,m) = Z;
        % Segmentation smoothing
        s2=0;
        while abs(s2 - smu2(k-kfrom+1,m))>eps*smu2(k-kfrom+1,m)
            Ms = arr2mat(Z',k);
            % step 5a in segmentation smoothing algorithm
            Nbkt = cell2mat(cellfun(@(x)conv(x,ones(2*b+1,1),'same'),mat2cell(Ms,ones(1,size(Ms,1)),[size(Ms,2)]),'UniformOutput',false))-Ms;
            % step 5b in segmentation smoothing algorithm
            [~,Z] = min((repmat(diag(Y'*Y),1,k)'-(W'*Y).^2)./(2*smu2(k-kfrom+1,m)*(J-1))-lambda*Nbkt,[],1);
            smu2(k-kfrom+1,m) = sum(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(T*(J-1));
            % step 7
            s2 = smu2(k-kfrom+1,m);
        end
    end
    c = c+1;
    mcv = min(smu2,[],2)'.*((J-1)^(-1)*(J-1-[kfrom:kto])).^(-2);
    [~,bestK] = min(mcv);
    if smu2(k-kfrom+1,m) < best_energy_so_far
        bestK = k;
        bestZ = Z;
        bestW = W;
        best_energy_so_far = smu2(k-kfrom+1,m);
        winner = m;
    end
end
k = bestK;
W = bestW;
Z = bestZ;

A = diag(Y'*W(:,Z));
sD2 = sum(diag(Y'*Y))/(T*(J-1));
R2 = 1 - best_energy_so_far/sD2;
W = bestW;
K = bestK;
end