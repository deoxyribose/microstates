function M=arr2mat(a,K)
% convert assignment array to matrix
% a is the assignmen 1*N
% M is matrix K*N
% K is the number of clusters if nargin==1 K=max(a)
N=length(a);
if nargin==1, K=max(a);end
M=zeros(K*N,1);
M(a+(0:(N-1))*K)=1;
M=reshape(M,K,N);