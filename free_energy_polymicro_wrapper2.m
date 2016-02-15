function [F] = free_energy_polymicro_wrapper2(gammas)
params = load('passedparams.mat');
X0 = rand(params.K,size(params.Y_train,2));
[W_train,X,M,beta,gamma1,gamma2,nits] = polymicro_smooth(params.Y_train,params.K,params.draw,gammas(1),gammas(2));
%nits
W = params.Sigma12invSigma22*W_train;
Wsq = W'*W;
diagww = diag(Wsq);
WY = W'*params.Y_test;
%mse = mean(mean((W_pred-params.w_test).^2));
%outparams = {gamma1,gamma2,beta,W,X,M};
[J,T] = size(params.Y_test);
logistic = @(x)1/(exp(-x)+1);
for i=1:3
    yrec = W*(M.*X);
    M = 0.95*M + 0.05*arrayfun(logistic,beta*X.*(-0.5*bsxfun(@times,X.*(1-2*M),diag(W'*W)) + W'*(params.Y_test-yrec)) + gamma1 + gamma2*([M(:,2:end) M(:,1)] + [M(:,end) M(:,1:end-1)]));
    
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
end
F1 = .5*J*T*log(beta/(2*pi))-beta/2*sum(sum((params.Y_test-W*(X.*M)).^2))-beta/2*sum(sum(diag(diag(W'*W))*(X.^2.*(M-M.^2))));
%F2 = log(gamma00) + sum(sum(M*log(gamma01/gamma00) + [M(:,end) M(:,1:end-1)]*log(gamma01/gamma00) + M.*[M(:,end) M(:,1:end-1)]*log((gamma00*gamma11)/(gamma01*gamma10))));
F3 = - sum(sum(M.*log(M) + (1-M).*log(1-M)));
F = double(F1+F3);
end