function K = fastrbf(x1,sig,varargin)

%function K = rbf(x1,sig,x2)
%
% Computes an rbf kernel matrix from the input coordinates
%
%INPUTS
% x1 =  a matrix containing all samples as rows
% sig = sigma, the kernel width; squared distances are divided by
%       squared sig in the exponent
% x2 =  (optional) a matrix containing all samples as rows, to make a
% rectangular kernel
%
%OUTPUTS
% K = the rbf kernel matrix ( = exp(-1/(2*sigma^2)*(coord*coord')^2) )
%
%
% For more info, see www.kernel-methods.net

switch nargin
    case 2
        n=size(x1,1);
        K=repmat(diag(x1*x1')',n,1)+repmat(diag(x1*x1'),1,n)-2*(x1*x1');
    case 3
        x2=varargin{1};
        n=size(x1,1);
        n2=size(x2,1);
        K=repmat(diag(x2*x2')',n,1)+repmat(diag(x1*x1'),1,n2)-2*(x1*x2');
    otherwise
        disp 'error, wrong number of arguments to rbf function'
end
K=exp((-1/(2*sig^2))*K);

%Author: Tijl De Bie, february 2003. Adapted: october 2004 (for speedup), 
%       Adapted again: Martin Kolar, on May 2012 (for unsquare kernels).
