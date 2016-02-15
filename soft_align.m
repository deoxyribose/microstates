function [est_perm,signs]=soft_align(h1,h2,H1,H2)
% This function estimates the alignment permutation matrix between
% the two sets of hidden unit vector sets H1, H2 of *same* size. And
% the relative signs to best match H1,H2
%
% Finds the sequence of best matches according to correlation coefficient
%
% INPUT
%       H1 master units K*N  (K number of hidden units, N sample size)
%       H2 array to be aligned with H1
% OUTPUT
%       est_perm  K-array with permutations
%       signs     K array with the sign factors (+/- 1)
%
% (c) Lars Kai Hansen July 2005
% edited by Franciszek Zdyb February 2016
%  IMM DTU
%


[K,N]=size(H1);
[K2,N2]=size(H2);
est_perm=randi(K,K,1);
one_to_K=zeros(1,K);
if K2~=K, disp('Error: Not same hidden dimension'), end
if N2~=N, disp('Error: Not same sample size'), end
COVM = corr(H2',H1');
% COVW = corr(h2,h1);
% 
% % find permutation based on h2, test on H2
% while any(one_to_K ~= 1:K)
%     [~,est_perm1] = max(abs(COVW));
%     [~,est_perm2] = max(abs(COVW),[],2);
%     [~,est_perm3] = max(abs(COVM));
%     [~,est_perm4] = max(abs(COVM),[],2);
%     COVM = COVM(est_perm,:);
%     [~,one_to_K] = max(abs(COVM));
%     signs = sign(diag(COVW(est_perm,:)));
% end
% %whos
for k=1:K,
    [dum,ind]=max(reshape(abs(COVM),K*K,1));
    ind1=rem(ind-1,K)+1;
    ind2=ceil(ind/K);
    % COVM(ind1,ind2) is max value
    est_perm(ind1)=ind2;
    signs(ind1)=sign(COVM(ind1,ind2));
    COVM(ind1,:)=zeros(1,K);
    COVM(:,ind2)=zeros(K,1);
end,
if numel(unique(est_perm)) ~= K
    missing = setdiff([1:K],est_perm);
    if numel(missing)>1
        est_perm = [1:K]';
        signs = ones(1,K);
    else
    [~,duplicate] = max(hist(est_perm,unique(est_perm)))
    est_perm(duplicate) = missing;
    COVM = corr(H2',H1');
    signs = sign(diag(COVM(est_perm,:)));
    end
end