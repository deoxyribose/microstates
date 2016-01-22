function [W,A,Z,K,clustering_algorithm] = basicNmicrostates(Y,kfrom,kto,MULTI,b,lambda)
clustering_algorithm = 'N-microstates';
[J,T] = size(Y);
for k=kfrom:kto
    best_energy_so_far = Inf;
    for m=1:MULTI
        disp(['Running N-microstates with ', num2str(k),' clusters for the ', num2str(m), 'th time'])
        % The basic N-microstate algorithm
        s2 = 0;
        eps = 10^(-10);
        % init
        W = datasample(Y,k,2); % take random frames
        W = bsxfun(@rdivide,W,sqrt(sum(W.^2,1)));
        % step 3
        [~,Z] = max((Y'*W).^2,[],2);
        c = 0;
        smu2 = sum(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(T*(J-1));
        while abs(s2 - smu2)>eps*smu2 % minimize orthogonal squared distance between each observed vector and corresponding microstates
            s2 = smu2;
            [~,Z] = max((Y'*W).^2,[],2);
            W = zeros(J,k);
            S = zeros(J,J,k);
            for i = 1:k;
                S(:,:,i) = Y(:,Z==i)*Y(:,Z==i)'; % sum covariances of all spatial maps for microstate i
                %W(:,i) = eig(Y(:,Z==i)*Y(:,Z==i)');
%                 for t=1:T
%                     if Z(t) == i
%                         S = S + Y(:,t)*Y(:,t)';
%                     end
%                 end
                [W(:,i),~] = eigs(S(:,:,i),1); % In what direction is the covariance largest? Eigs returns normalized eigenvector corresp. to largest eigenvalue
            end
            % W = bsxfun(@rdivide,W,sqrt(sum(W.^2,1)));
            smu2 = sum(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(T*(J-1));
            c = c+1;
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
        if smu2 < best_energy_so_far
            bestK = k;
            bestZ = Z;
            bestW = W;
            best_energy_so_far = smu2
        end
        sD2 = sum(diag(Y'*Y))/(T*(J-1));
        R2 = 1 - smu2/sD2;
    end
end
k = bestK;
W = bestW;
Z = bestZ;
% segmentation smoothing
s2=0;
% while abs(s2 - smu2)>eps*smu2
%     for t=(1+b):(T-b)
%         for i=1:k
%             Nbkt = sum(Z((t-b):(t+b)) == i);
%         end
%         [~,Z] = min(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(2*smu2*(J-1)-lambda*Nbkt);
%     end
%     smu2 = sum(diag(Y'*Y)-diag((W(:,Z)'*Y).^2))/(T*(J-1));
% end
% a = zeros(k,J);
% for t=1:T;
%     kappa = Z(t);
%     for i=1:k
%         if kappa == i
%             a(kappa,t) = Y(:,t)'*W(:,kappa);
%         end
%     end
% end
A = diag(Y'*W(:,Z));
sD2 = sum(diag(Y'*Y))/(T*(J-1));
smu2;
R2 = 1 - smu2/sD2;
Z = Z';
W = bestW;
K = bestK;
end