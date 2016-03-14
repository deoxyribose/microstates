% Copyright (C) 2015  Franciszek Zdyb
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [OUTEEG, com] = pop_getmicrostates( INEEG, varargin);

% SUMMARY:
%    pop_getmicrostates clusters a <subset> of <INEEG.data> into a varying number of
%    microstates, from <Kfrom> to <Kto>, using
%    <clustering_algorithm>, and selects the best number according to the
%    Krzanovski-Lai criterion.

% INPUTS:
%   <subset> is one of {all data, downsampled data, local maxima in global field power}
%   <clustering_algorithm> is one of {FastICA, k-means, agglomerative clustering}

% OUTPUTS:
com = ''; % this initialization ensure that the function will return something
% if the user press the cancel button
OUTEEG = INEEG;

% varargin = [subset, Kfrom, Kto, clustering_algorithm, draw,
% lambda, chronos, W]
if size(size(OUTEEG.data),2) > 2
    disp(['Data has more dimensions than nbchan x pnts, the first 3 dims are : ', num2str(size(OUTEEG.data,1)),' x ',num2str(size(OUTEEG.data,2)),' x ',num2str(size(OUTEEG.data,3))])
    disp('Do something like OUTEEG.data = OUTEEG.data(:,:,i_epoch) and press Continue')
    disp('or!')
    disp('...')
    disp('You could do something like this: [OUTEEG1, com] = pop_getmicrostates( INEEG, varargin); OUTEEG.data = OUTEEG.data(:,:,i_epoch1); [OUTEEG21, com] = pop_getmicrostates( INEEG, varargin); OUTEEG.data = OUTEEG.data(:,:,i_epoch2);')
    disp('...and press continue')
    keyboard
end
% the defaults
OUTEEG.subset = 1;
OUTEEG.Kfrom = 3;
OUTEEG.Kto = 7;
OUTEEG.clustering_algorithm = 6;
startframe = 1;
endframe = OUTEEG.pnts;
chronos = 1;
draw = 0;
MULTI = 10;
W_init = 0;
for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~isstr(Param)
        error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'clustering_algorithm'
            OUTEEG.clustering_algorithm = Value;
        case 'subset'
            OUTEEG.subset = Value;
        case 'kfrom'
            OUTEEG.Kfrom = Value;
        case 'kto'
            OUTEEG.Kto = Value;
        case 'draw'
            draw = Value;
        case 'multi'
            MULTI = Value;
        case 'chronos'
            chronos = Value;
        case 'verbose'
            verbose = Value;
        case 'w_init'
            OUTEEG.W_init = Value;
        otherwise
            error(['Unknown input parameter ''' Param ''' ???'])
    end
end

if nargin < 3
    [~, ~, ~, structout] = inputgui( 'geometry', [2 3 2 3 1], ...
        'geomvert', [], 'uilist', { ...
        { 'style', 'text', 'string', [ 'Choose subset of data to run clustering on:' 10 10 10] }, ...
        { 'style', 'listbox', 'string', 'All data|Peak topographies (requires that pop_get_clean_topographies has been run)|Downsampled data' 'tag' 'subset' }, ...
        { 'style', 'text', 'string', [ 'Choose range of clusters to test:' 10 ] }, ...
        { 'style', 'edit', 'string', '4' 'tag' 'Kfrom'} , ...
        { 'style', 'edit', 'string', '20' 'tag' 'Kto'} , ...
        { 'style', 'text', 'string', [ 'Choose clustering algorithm:' 10 10 10 10] }, ...
        { 'style', 'listbox', 'string', 'ICA|Kmeans|Agglomerative|N-microstate|DPmeans|Variational Microstates|Polymicrostates' 'tag' 'algorithm' }, ...
        { 'style', 'text', 'string', [ 'Choose subset of data (default is all)' 10 ] }, ...
        { 'style', 'edit', 'string', '1' 'tag' 'start'} , ...
        { 'style', 'edit', 'string', num2str(OUTEEG.pnts) 'tag' 'end'} , ...
        { 'style', 'checkbox', 'string', 'Use chronological order for clustering' 'tag' 'chronos', 'Value',1 }, ...
        },'title','Microstate analysis' );
    OUTEEG.subset = structout.subset;
    OUTEEG.Kfrom = str2num(structout.Kfrom);
    OUTEEG.Kto = str2num(structout.Kto);
    OUTEEG.clustering_algorithm = structout.algorithm;
    % User can further constrain the data between frame nr. <start> and
    % <end>
    startframe = str2num(structout.start);
    endframe = str2num(structout.end);
    % User can toggle temporal smoothing
    chronos = structout.chronos;
    draw = 1;
end;

% Now process inputs that aren't algorithm-specific.
switch OUTEEG.subset
    case 1
        OUTEEG.peakidx = 1:size(OUTEEG.data,2);
        subset = 'all';
    case 2
        subset = ['GFP-peak topographies'];
        if ~isfield(OUTEEG, 'peakidx')
            OUTEEG = pop_get_clean_topographies(OUTEEG, 1); %get GFP peaks
        end
    case 3
        promptstr    = { 'Keep every nth frame (default is every 2nd)' };
        inistr       = { '2' };
        result       = inputdlg( promptstr, 'Title of window', 1,  inistr);
        if isempty(result) return; end;
        everynth   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
        OUTEEG.peakidx = 1:everynth:size(OUTEEG.data,2);
        subset = [num2str(length(OUTEEG.peakidx)/OUTEEG.pnts*100), ' %'];
end
rbfsigma = 3287.86; alpha = 3748.3; % found with bayesopthyperparams.m, but probably suboptimal
OUTEEG.peakidx = OUTEEG.peakidx(OUTEEG.peakidx>=startframe);
OUTEEG.peakidx = OUTEEG.peakidx(OUTEEG.peakidx<=endframe);
Y = OUTEEG.data(:,OUTEEG.peakidx);
if chronos
    Y = [Y; (1:size(Y,2))/alpha]; % temporal smoothing is done by adding a scaled timeline, [1,2,3...,OUTEEG.pnts], as a feature
end
% standardize data, i.e. substract mean and divide by standard deviation
OUTEEG.data_std = std(Y,[],2);
OUTEEG.data_mean = mean(Y,2);
Y=bsxfun(@times,bsxfun(@minus,Y,OUTEEG.data_mean),1./OUTEEG.data_std);


number_of_ks = OUTEEG.Kto-OUTEEG.Kfrom+1;
label = zeros(size(Y,2),number_of_ks);
N = length(OUTEEG.peakidx);
D = OUTEEG.nbchan;
krange = OUTEEG.Kfrom:OUTEEG.Kto;
disp(['Restarts ', num2str(MULTI), ' times']);
% Now go to chosen algorithm
switch OUTEEG.clustering_algorithm
    case 1 % ICA
        clustering_algorithm = 'ICA';
        if OUTEEG.Kto~=OUTEEG.Kfrom
            for k=krange
                [~, W] = fastica(Y,'numOfIC',k);
                S = W*Y;
                % calculate microstate sequence, as per Yuan 2012
                [~, lbl] = max(abs(S'),[],2);
                label(:,k-OUTEEG.Kfrom+1) = lbl';
                disp('Done with ICA')
            end
            % KL
            [KL, w] = getKL(Y, label,OUTEEG.Kfrom);
            if draw
                figure('name','KL-criterion')
                subplot(211)
                plot(krange,w)
                ylabel('Dispersion')
                subplot(212)
                plot(krange(2:end-1),KL)
                ylabel('Krzanovski-Lai criterion')
            end
            [~,best_k] = max(KL); % argmax
            
            OUTEEG.K = best_k + OUTEEG.Kfrom;
        else
            OUTEEG.K = OUTEEG.Kto;
        end
        [A, W] = fastica(Y,'numOfIC',OUTEEG.K);
        S = W*Y;
        
        vars = get_component_energies(A,S);
        [OUTEEG.vars, order] = sort(vars);
        OUTEEG.W = A(:,order);
        OUTEEG.S = S(order,:);
        
        [~, OUTEEG.Z] = max(abs(S'),[],2);
        OUTEEG.Z = OUTEEG.Z';
        
        % global map dissimilarity
        if chronos
            OUTEEG.meanGMD = mean(GMD(Y(1:(size(Y,1)),:),A(:,OUTEEG.Z),size(Y,1)));
        else
            OUTEEG.meanGMD = mean(GMD(Y,A(:,OUTEEG.Z),OUTEEG.nbchan));
        end
        if draw
            if size(Y,2) <= 5000
                show_clusters(Y,OUTEEG.Z,OUTEEG.meanGMD,OUTEEG.K);
            else
                show_clusters(Y(:,1:5000),OUTEEG.Z(1:5000),OUTEEG.meanGMD,OUTEEG.K);
            end
        end
    case 2 % kmeans
        clustering_algorithm = 'kmeans';
        label = zeros(size(Y,2),number_of_ks);
        %K = double(fastrbf(Y',rbfsigma)); % double() because otherwise bug in knkmeans in line 22
        %K = doubles(kernelmatrix('poly',Y,Y,0.01,0,4));
        if OUTEEG.Kto~=OUTEEG.Kfrom
            for k=krange
                best_energy_so_far = Inf;
                c = 0;
                for m=1:MULTI
                    %[lbl, energy] = knkmeans(K,k);
                    [W,lbl,~,energy] = kmeans_fast(Y',k);
                    if energy < best_energy_so_far
                        bestW = W;
                        label(:,k-OUTEEG.Kfrom+1) = lbl;
                        best_energy_so_far = energy;
                        c = c+1;
                    end
                end
                %show_clusters(Y,label,energy,k);
                nclusters = numel(unique(label(:,k-OUTEEG.Kfrom+1)));
                if nclusters < k
                    disp(['Couldnt estimate more than ', num2str(nclusters), ' clusters'])
                    break
                end
            end
            % KL
            [KL, w] = getKL(Y, label,OUTEEG.Kfrom);
            if draw
                figure('name','KL-criterion')
                subplot(211)
                plot(krange,w)
                ylabel('Dispersion')
                subplot(212)
                plot(krange(2:end-1),KL)
                ylabel('Krzanovski-Lai criterion')
            end
            [~,best_k] = max(KL);
            OUTEEG.K = best_k + OUTEEG.Kfrom;
            OUTEEG.Z = label(:,best_k+1)';
        else
            OUTEEG.K = OUTEEG.Kto;
            c = 0;
            best_energy_so_far = Inf;
            for m=1:MULTI
                [W,lbl,~,energy] = kmeans_fast(Y',OUTEEG.K);
                if energy < best_energy_so_far
                    best_energy_so_far = energy;
                    c = c+1;
                end
            end
            OUTEEG.Z = lbl';
        end
        %         A = zeros(OUTEEG.nbchan,OUTEEG.K);
        %         for iMCRST=1:OUTEEG.K
        %             A(:,iMCRST) = mean(OUTEEG.data(:,OUTEEG.Z==iMCRST),2);
        %         end
        OUTEEG.W = W';
        OUTEEG.A = 0;
%         if chronos
%             OUTEEG.meanGMD = mean(GMD(Y(1:(size(Y,1)-1),:),OUTEEG.W(:,OUTEEG.Z),size(Y,1)-1));
%         else
%             OUTEEG.meanGMD = mean(GMD(Y,OUTEEG.W(:,OUTEEG.Z),OUTEEG.nbchan));
%         end
        if draw
            if size(Y,2) <= 5000
                show_clusters(Y,OUTEEG.Z,0,OUTEEG.K);
            end
        end
        disp('Done with kmeans')
    case 3
        clustering_algorithm = 'agglomerative clustering';
        disp('Running agglomerative clustering...')
        distmat = pdist(Y',@GMD, OUTEEG.nbchan);
        Z = PHA_Clustering(squareform(distmat));
        
        if draw
            if size(Y,2) < 2000
                figure('name','dendrogram')
                dendrogram(Z);
            end
        end
        if OUTEEG.Kto~=OUTEEG.Kfrom
            for k=krange
                label(:,k-OUTEEG.Kfrom+1) = cluster(Z,'maxclust',k);
            end
            % KL
            [KL, w] = getKL(Y, label,OUTEEG.Kfrom);
            if draw
                figure('name','KL-criterion')
                subplot(211)
                plot(krange,w)
                ylabel('Dispersion')
                subplot(212)
                plot(krange(2:end-1),KL)
                ylabel('Krzanovski-Lai criterion')
            end
            [~,best_k] = max(KL);
            OUTEEG.Z = label(:,best_k+1)';
            OUTEEG.K = best_k + OUTEEG.Kfrom;
        else
            OUTEEG.K = OUTEEG.Kto;
            OUTEEG.Z = cluster(Z,'maxclust',OUTEEG.K)';
        end
        W = zeros(OUTEEG.nbchan,OUTEEG.K);
        for iMCRST=1:OUTEEG.K
            W(:,iMCRST) = mean(OUTEEG.data(:,OUTEEG.Z==iMCRST),2);
        end
        OUTEEG.W = W;
        OUTEEG.A = 0;
%         if chronos
%             OUTEEG.meanGMD = mean(GMD(Y(1:(size(Y,1)-1),:),A(:,OUTEEG.Z),size(Y,1)-1));
%         else
%             OUTEEG.meanGMD = mean(GMD(Y,A(:,OUTEEG.Z),OUTEEG.nbchan));
%         end
        if draw
            if size(Y,2) <= 5000
                show_clusters(Y,OUTEEG.Z,0,OUTEEG.K);
            end
        end
        disp('Done with agglomerative clustering')
    case 4
        b = 3;
        lambda = 5;
        [OUTEEG.W,OUTEEG.A,OUTEEG.Z,~,~,~,~,OUTEEG.K,clustering_algorithm] = basicNmicrostates(double(Y),OUTEEG.Kfrom,OUTEEG.Kto,MULTI,b,lambda,10^(-6),verbose);
        OUTEEG.A = repmat(OUTEEG.A',OUTEEG.K,1);
    case 5
        clustering_algorithm = 'DP-means';
        Y = Y(1:end-1,:);
        %lambda = 257; % found with bayesopt
        lambda = 150;
        k = 1;
        z = ones(length(Y),1);
        oldz = z;
        mu = mean(Y,2);
        converged = 0;
        oldobj = Inf;
        while ~converged
            d = zeros(length(Y),k);
            for i=1:length(Y)
                for c=1:k
                    dic = Y(:,i)-mu(:,c);
                    d(i,c) = dic'*dic;
                    if isnan(d(i,c))
                        disp('lol');
                    end
                end
                if min(d(i,:))>lambda
                    k = k+1;
                    z(i) = k;
                    mu = [mu Y(:,i)];
                else
                    [~, z(i)] = min(d(i,:));
                end
            end
            if length(unique(z)) < max(z)
                z
                z = z-1;
                k = k-1;
                mu = mu(:,2:end);
            end
            for c=1:k
                mu(:,c) = mean(Y(:,z==c),2);
                if isnan(mu(:,c))
                    disp('lol');
                end
            end
            if oldz == z
                converged = 1;
            else
                oldz = z;
            end
            %             obj = sum(sum(d)) + lambda*k;
            %             rel_improvement = abs(1 - oldobj/obj);
            %             converged = rel_improvement < 0.0001; % change to no changes in assignment
            %             oldobj = obj;
        end
        OUTEEG.K = k;
        OUTEEG.W = mu;
        OUTEEG.Z = z';
        disp('Done with DP-means')
    case 6
        clustering_algorithm = 'Variational Microstates';
        if chronos
            Y = Y(1:end-1,:);
        end
        if OUTEEG.Kto~=OUTEEG.Kfrom
            disp(['Model selection not implemented yet, proceeding with ', num2str(OUTEEG.Kfrom)])
            OUTEEG.K = OUTEEG.Kfrom;
        end
        draw = 0;
        K = OUTEEG.Kfrom;
        alpha = 2.5;
        gamma2 = 5;
        G = 0; %trigger
        max_nits = 5000;
        kfolds = 1;
        learn_rate_init = 0.01;
        learn_decay = 0.0001;
        W_0 = W_init;
        verbose = 1;
        [W_var,Mu,M,sigma2,allZs,beta_var,free_energy,recon_error,m_winner,nits,varargin] = variational_microstates_smooth(Y,K,draw,alpha,gamma2,MULTI,max_nits,G,learn_rate_init,learn_decay,verbose);
        disp(['Converged in ', num2str(nits), ' iterations'])
        [~,OUTEEG.Z] = max(M,[],1);
        OUTEEG.K = K;
        OUTEEG.W = W_var;
        OUTEEG.A = Mu;
        OUTEEG.energy = free_energy;
    case 7
        clustering_algorithm = 'Polymicrostates';
        if chronos
            Y = Y(1:end-1,:);
        end
        if OUTEEG.Kto~=OUTEEG.Kfrom
            OUTEEG.K = size(OUTEEG.Y,1); % J
            disp(['Looking for K = nbchan microstates, K = ', num2str(OUTEEG.K)])
        else 
            OUTEEG.K = OUTEEG.Kto;
        end
        draw = 0;
        gamma1 = -50;
        gamma2 = 20;
        max_nits = 2000;
        learn_rate_init = 0.0005;
        learn_rate_decay = 0.00001;
        verbose = 1;
        [W,X,M,~,~,~,~,~,~,~,~] = polymicro_smooth(Y,'K',OUTEEG.K,'draw',draw,'MULTI',MULTI,'max_nits',max_nits,'learn_rate_init',learn_rate_init,'learn_rate_decay',learn_rate_decay,'gamma1',gamma1,'gamma2',gamma2,'verbose',verbose);
        OUTEEG.Z = M>0.5;
        OUTEEG.W = W;
        OUTEEG.A = X;
    otherwise
        disp('No.')
end
disp(['Estimated ', num2str(OUTEEG.K), ' microstates, using ', clustering_algorithm, ' on ', subset, ' of the data']);

% return the string command
% -------------------------
com = sprintf('pop_getmicrostates( %s, %d, [%s] );', inputname(1), int2str(OUTEEG.subset), int2str(OUTEEG.Kfrom), int2str(OUTEEG.Kto), int2str(OUTEEG.clustering_algorithm));

return;
