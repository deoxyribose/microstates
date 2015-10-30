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

function [OUTEEG, com] = pop_getmicrostates( INEEG, subset, nMCRSTSfrom, nMCRSTSto, clustering_algorithm, draw);

% SUMMARY:
% pop_getmicrostates clusters a <subset> of <INEEG.data> into a varying number of
% microstates, from <nMCRSTSfrom> to <nMCRSTSto>, using
% <clustering_algorithm>, and selects the best number according to the
% Krzanovski-Lai criterion.

% <subset> is one of {all data, downsampled data, local maxima in global field power}
% <clustering_algorithm> is one of {FastICA, k-means, agglomerative
% clustering}

com = ''; % this initialization ensure that the function will return something
% if the user press the cancel button
OUTEEG = INEEG;

% GUI asking for analysis arguments
% -------------
if nargin < 3
    [~, ~, ~, structout] = inputgui( 'geometry', [2 3 2 3 1], ...
        'geomvert', [], 'uilist', { ...
        { 'style', 'text', 'string', [ 'Choose subset of data to run clustering on:' 10 10 10] }, ...
        { 'style', 'listbox', 'string', 'All data|Peak topographies (requires that pop_get_clean_topographies has been run)|Downsampled data' 'tag' 'subset' }, ...
        { 'style', 'text', 'string', [ 'Choose range of clusters to test:' 10 ] }, ...
        { 'style', 'edit', 'string', '4' 'tag' 'nMCRSTSfrom'} , ...
        { 'style', 'edit', 'string', '20' 'tag' 'nMCRSTSto'} , ...
        { 'style', 'text', 'string', [ 'Choose clustering algorithm:' 10 10 10 10] }, ...
        { 'style', 'listbox', 'string', 'ICA|Kmeans|Agglomerative|N-microstate|DPmeans' 'tag' 'algorithm' }, ...
        { 'style', 'text', 'string', [ 'Choose subset of data (default is all)' 10 ] }, ...
        { 'style', 'edit', 'string', '1' 'tag' 'start'} , ...
        { 'style', 'edit', 'string', num2str(OUTEEG.pnts) 'tag' 'end'} , ...
        { 'style', 'checkbox', 'string', 'Use chronological order for clustering' 'tag' 'chronos', 'Value',1 }, ...
        },'title','Microstate analysis' );
    OUTEEG.subset = structout.subset;
    OUTEEG.nMCRSTSfrom = str2num(structout.nMCRSTSfrom);
    OUTEEG.nMCRSTSto = str2num(structout.nMCRSTSto);
    OUTEEG.clustering_algorithm = structout.algorithm;
    % User can further constrain the data between frame nr. <start> and
    % <end>
    startframe = str2num(structout.start);
    endframe = str2num(structout.end);
    % User can toggle temporal smoothing
    chronos = structout.chronos;
    draw = 1;
else
    % If function is run without GUI and no params are provided, these are
    % the defaults
    OUTEEG.subset = subset;
    OUTEEG.nMCRSTSto = nMCRSTSto;
    OUTEEG.nMCRSTSfrom = nMCRSTSfrom;
    OUTEEG.clustering_algorithm = clustering_algorithm;
    startframe = 1;
    endframe = OUTEEG.pnts;
    chronos = 1;
    draw = 0;
end;
% How many times to restart k-means
MULTI = 100;
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
% Second GUI asks for MULTI
% if OUTEEG.clustering_algorithm == 2
%     promptstr    = { 'How many restarts of kmeans to run:' };
%     inistr       = { '7' };
%     result       = inputdlg( promptstr, 'Title of window', 1,  inistr);
%     if length( result ) == 0 return; end;
%
%     MULTI   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
% end
% ---------------------------------------------------
rbfsigma = 3287.86; alpha = 3748.3; % found with bayesopthyperparams.m, but probably suboptimal
OUTEEG.peakidx = OUTEEG.peakidx(OUTEEG.peakidx>=startframe);
OUTEEG.peakidx = OUTEEG.peakidx(OUTEEG.peakidx<=endframe);
signal = OUTEEG.data(:,OUTEEG.peakidx);
if chronos
    signal = [signal; (1:size(signal,2))/alpha]; % temporal smoothing is done by adding a scaled timeline, [1,2,3...,OUTEEG.pnts], as a feature
end
% standardize data, i.e. substract mean and divide by standard deviation
signal=bsxfun(@rdivide,bsxfun(@minus,signal,mean(signal,2)),std(signal,[],2));
%     if PCAs
%         [U,Sigma,V]=svdecon(X');
%         signal_pca = U*Sigma;
%         signal = signal_pca';
%         %S2 = Sigma.*Sigma;
%         %find(cumsum(diag(S2))/trace(S2)<0.99);
%     else
%         signal = X;
%     end
number_of_ks = OUTEEG.nMCRSTSto-OUTEEG.nMCRSTSfrom+1;
label = zeros(size(signal,2),number_of_ks);
N = length(OUTEEG.peakidx);
D = OUTEEG.nbchan;
krange = OUTEEG.nMCRSTSfrom:OUTEEG.nMCRSTSto;
% GUI
switch OUTEEG.clustering_algorithm
    case 1 % ICA
        clustering_algorithm = 'ICA';
        if OUTEEG.nMCRSTSto~=OUTEEG.nMCRSTSfrom
            for k=krange
                [~, W] = fastica(signal,'numOfIC',k);
                S = W*signal;
                % calculate microstate sequence, as per Yuan 2012
                [~, lbl] = max(abs(S'),[],2);
                label(:,k-OUTEEG.nMCRSTSfrom+1) = lbl';
                disp('Done with ICA')
            end
            % KL
            [KL, w] = getKL(signal, label,OUTEEG.nMCRSTSfrom);
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
            
            OUTEEG.nMCRSTS = best_k + OUTEEG.nMCRSTSfrom;
        else
            OUTEEG.nMCRSTS = OUTEEG.nMCRSTSto;
        end
        %         OUTEEG.idx = label(:,best_k+1)';
        %         A = zeros(OUTEEG.nbchan,OUTEEG.nMCRSTS);
        %         for iMCRST=1:OUTEEG.nMCRSTS
        %             A(:,iMCRST) = mean(OUTEEG.data(:,OUTEEG.idx==iMCRST),2);
        %         end
        %
        %         OUTEEG.A = A;
        [A, W] = fastica(signal,'numOfIC',OUTEEG.nMCRSTS);
        S = W*signal;
        
        vars = get_component_energies(A,S);
        [OUTEEG.vars, order] = sort(vars);
        OUTEEG.A = A(:,order);
        OUTEEG.S = S(order,:);
        
        [~, OUTEEG.idx] = max(abs(S'),[],2);
        OUTEEG.idx = OUTEEG.idx';
        
        % global map dissimilarity
        if chronos
            OUTEEG.meanGMD = mean(GMD(signal(1:(size(signal,1)),:),A(:,OUTEEG.idx),size(signal,1)));
        else
            OUTEEG.meanGMD = mean(GMD(signal,A(:,OUTEEG.idx),OUTEEG.nbchan));
        end
        if draw
            if size(signal,2) <= 5000
                show_clusters(signal,OUTEEG.idx,OUTEEG.meanGMD,OUTEEG.nMCRSTS);
            else
                show_clusters(signal(:,1:5000),OUTEEG.idx(1:5000),OUTEEG.meanGMD,OUTEEG.nMCRSTS);
            end
        end
    case 2 % kmeans
        clustering_algorithm = 'kmeans';
        label = zeros(size(signal,2),number_of_ks);
        K = double(fastrbf(signal',rbfsigma)); % double() because otherwise bug in knkmeans in line 22
        %K = doubles(kernelmatrix('poly',signal,signal,0.01,0,4));
        if OUTEEG.nMCRSTSto~=OUTEEG.nMCRSTSfrom
            for k=krange
                best_energy_so_far = Inf;
                c = 0;
                for m=1:MULTI
                    [lbl, energy] = knkmeans(K,k);
                    if energy < best_energy_so_far
                        label(:,k-OUTEEG.nMCRSTSfrom+1) = lbl;
                        best_energy_so_far = energy;
                        c = c+1;
                    end
                end
                %show_clusters(signal,label,energy,k);
                nclusters = numel(unique(label(:,k-OUTEEG.nMCRSTSfrom+1)));
                if nclusters < k
                    disp(['Couldnt estimate more than ', num2str(nclusters), ' clusters'])
                    break
                end
            end
            % KL
            [KL, w] = getKL(signal, label,OUTEEG.nMCRSTSfrom);
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
            OUTEEG.nMCRSTS = best_k + OUTEEG.nMCRSTSfrom;
            OUTEEG.idx = label(:,best_k+1)';
        else
            OUTEEG.nMCRSTS = OUTEEG.nMCRSTSto;
            c = 0;
            best_energy_so_far = Inf;
            for m=1:MULTI
                [lbl, energy] = knkmeans(K,OUTEEG.nMCRSTS);
                if energy < best_energy_so_far
                    best_energy_so_far = energy;
                    c = c+1;
                end
            end
            OUTEEG.idx = lbl;
        end
        A = zeros(OUTEEG.nbchan,OUTEEG.nMCRSTS);
        for iMCRST=1:OUTEEG.nMCRSTS
            A(:,iMCRST) = mean(OUTEEG.data(:,OUTEEG.idx==iMCRST),2);
        end
        OUTEEG.A = A;
        if chronos
            OUTEEG.meanGMD = mean(GMD(signal(1:(size(signal,1)-1),:),A(:,OUTEEG.idx),size(signal,1)-1));
        else
            OUTEEG.meanGMD = mean(GMD(signal,A(:,OUTEEG.idx),OUTEEG.nbchan));
        end
        if draw
            if size(signal,2) <= 5000
                show_clusters(signal,OUTEEG.idx,0,OUTEEG.nMCRSTS);
            end
        end
        disp('Done with kmeans')
    case 3
        clustering_algorithm = 'agglomerative clustering';
        disp('Running agglomerative clustering...')
        distmat = pdist(signal',@GMD, OUTEEG.nbchan);
        Z = PHA_Clustering(squareform(distmat));
        
        if draw
            if size(signal,2) < 2000
                figure('name','dendrogram')
                dendrogram(Z);
            end
        end
        if OUTEEG.nMCRSTSto~=OUTEEG.nMCRSTSfrom
            for k=krange
                label(:,k-OUTEEG.nMCRSTSfrom+1) = cluster(Z,'maxclust',k);
            end
            % KL
            [KL, w] = getKL(signal, label,OUTEEG.nMCRSTSfrom);
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
            OUTEEG.idx = label(:,best_k+1)';
            OUTEEG.nMCRSTS = best_k + OUTEEG.nMCRSTSfrom;
        else
            OUTEEG.nMCRSTS = OUTEEG.nMCRSTSto;
            OUTEEG.idx = cluster(Z,'maxclust',OUTEEG.nMCRSTS);
        end
        A = zeros(OUTEEG.nbchan,OUTEEG.nMCRSTS);
        for iMCRST=1:OUTEEG.nMCRSTS
            A(:,iMCRST) = mean(OUTEEG.data(:,OUTEEG.idx==iMCRST),2);
        end
        OUTEEG.A = A;
        if chronos
            OUTEEG.meanGMD = mean(GMD(signal(1:(size(signal,1)-1),:),A(:,OUTEEG.idx),size(signal,1)-1));
        else
            OUTEEG.meanGMD = mean(GMD(signal,A(:,OUTEEG.idx),OUTEEG.nbchan));
        end
        if draw
            if size(signal,2) <= 5000
                show_clusters(signal,OUTEEG.idx,0,OUTEEG.nMCRSTS);
            end
        end
        disp('Done with agglomerative clustering')
    case 4
        clustering_algorithm = 'N-microstates';
        signal = signal(1:end-1,:);
        for k=OUTEEG.nMCRSTSfrom:OUTEEG.nMCRSTSto
            k
            best_energy_so_far = Inf;
            for m=1:MULTI
                % The basic N-microstate algorithm
                s2 = 0;
                eps = 10^(-10);
                % init
                Gamma = datasample(signal,k,2); % take random frames
                Gamma = bsxfun(@rdivide,Gamma,sqrt(sum(Gamma.^2,1)));
                % step 3
                [~,OUTEEG.idx] = max((signal'*Gamma).^2,[],2);
                c = 0;
                smu2 = 0;
                while abs(s2 - smu2)>eps*smu2
                    s2 = smu2;
                    [~,OUTEEG.idx] = max((signal'*Gamma).^2,[],2);
                    Gamma = zeros(OUTEEG.nbchan,k);
                    S = zeros(OUTEEG.nbchan,OUTEEG.nbchan);
                    for i = 1:k;
                        %S(:,:,i) = signal(:,OUTEEG.idx==i)*signal(:,OUTEEG.idx==i)';
                        %Gamma(:,i) = eig(signal(:,OUTEEG.idx==i)*signal(:,OUTEEG.idx==i)');
                        for t=1:OUTEEG.pnts
                            if OUTEEG.idx(t) == i
                                S = S + signal(:,t)*signal(:,t)';
                            end
                        end
                        Gamma(:,i) = eig(S);
                    end
                    Gamma = bsxfun(@rdivide,Gamma,sqrt(sum(Gamma.^2,1)));
                    smu2 = sum(diag(signal'*signal)-diag((Gamma(:,OUTEEG.idx)'*signal).^2))/(OUTEEG.pnts*(OUTEEG.nbchan-1));
                    c = c+1
                end
                
                a = zeros(k,OUTEEG.nbchan);
                for t=1:OUTEEG.pnts;
                    kappa = OUTEEG.idx(t);
                    for i=1:k
                        if kappa == i
                            a(kappa,t) = signal(:,t)'*Gamma(:,kappa);
                        end
                    end
                end
                if smu2 < best_energy_so_far
                    OUTEEG.nMCRSTS = k;
                    OUTEEG.A = Gamma;
                    best_energy_so_far = smu2;
                end
                sD2 = sum(diag(signal'*signal))/(OUTEEG.pnts*(OUTEEG.nbchan-1));
                R2 = 1 - smu2/sD2;
            end
        end
        % segmentation smoothing
        b = 5;
        lambda = 7;
        Lambda = OUTEEG.idx;
        e = sum(diag(signal'*signal)-diag((Gamma(:,OUTEEG.idx)'*signal).^2))/(OUTEEG.pnts*(OUTEEG.nbchan-1));
        while abs(s2 - smu2)>eps*smu2
            for t=(1+b):(OUTEEG.pnts-b)
                for i=1:k
                    Nbkt = sum(L((t-b):(t+b)) == i);
                end
                [~,Lambda] = min(diag(signal'*signal)-diag((Gamma(:,OUTEEG.idx)'*signal).^2))/(2*e*(OUTEEG.nbchan-1)-lambda*Nbkt);
            end
            OUTEEG.idx = Lambda;
            smu2 = sum(diag(signal'*signal)-diag((Gamma(:,OUTEEG.idx)'*signal).^2))/(OUTEEG.pnts*(OUTEEG.nbchan-1));
        end
        a = zeros(k,OUTEEG.nbchan);
        for t=1:OUTEEG.pnts;
            kappa = OUTEEG.idx(t);
            for i=1:k
                if kappa == i
                    a(kappa,t) = signal(:,t)'*Gamma(:,kappa);
                end
            end
        end
        sD2 = sum(diag(signal'*signal))/(OUTEEG.pnts*(OUTEEG.nbchan-1));
        R2 = 1 - smu2/sD2;
        OUTEEG.idx = OUTEEG.idx';
    case 5
        clustering_algorithm = 'DP-means';
        signal = signal(1:end-1,:);
        lambda = 100;
        k = 1;
        z = ones(length(signal),1);
        mu = mean(signal,2);
        converged = 0;
        oldobj = Inf;
        while ~converged
            d = zeros(length(signal),k);
            for i=1:length(signal)
                for c=1:k
                    dic = signal(:,i)-mu(:,c);
                    d(i,c) = dic'*dic;
                    if isnan(d(i,c))
                        disp('lol');
                    end
                end
                if min(d(i,:))>lambda
                    k = k+1;
                    z(i) = k;
                    mu = [mu signal(:,i)];
                else
                    [~, z(i)] = min(d(i,:));
                end
            end
            if length(unique(z)) < max(z)
                z = z-1;
                k = k-1;
                mu = mu(:,2:end);
            end
            for c=1:k
                mu(:,c) = mean(signal(:,z==c),2);
                if isnan(mu(:,c))
                    disp('lol');
                end
            end
            obj = sum(sum(d)) + lambda*k;
            rel_improvement = abs(1 - oldobj/obj);
            converged = rel_improvement < 0.0001; % change to no changes in assignment
            oldobj = obj;
        end
        OUTEEG.nMCRSTS = k;
        OUTEEG.A = mu;
        OUTEEG.idx = z';
        disp('Done with DP-means')
    otherwise
        disp('No.')
end
disp(['Estimated ', num2str(OUTEEG.nMCRSTS), ' microstates, using ', clustering_algorithm, ' on ', subset, ' of the data']);

% return the string command
% -------------------------
com = sprintf('pop_getmicrostates( %s, %d, [%s] );', inputname(1), int2str(OUTEEG.subset), int2str(OUTEEG.nMCRSTSfrom), int2str(OUTEEG.nMCRSTSto), int2str(OUTEEG.clustering_algorithm));

return;
