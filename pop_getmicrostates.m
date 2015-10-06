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

function [OUTEEG, com] = pop_getmicrostates( INEEG, subset, nMCRSTSfrom, nMCRSTSto, clustering_algorithm);

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
        { 'style', 'listbox', 'string', 'ICA|Kmeans|Agglomerative' 'tag' 'algorithm' }, ...
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
end;
% How many times to restart k-means
MULTI = 1;
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
OUTEEG.peakidx = OUTEEG.peakidx(OUTEEG.peakidx>startframe);
OUTEEG.peakidx = OUTEEG.peakidx(OUTEEG.peakidx<endframe);
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
        for k=krange
            [A, W] = fastica(signal,'numOfIC',k);
            S = W*signal;
            % calculate microstate sequence, as per Yuan 2012
            [~, lbl] = max(abs(S'),[],2);
            label(:,k-OUTEEG.nMCRSTSfrom+1) = lbl';
            disp('Done with ICA')
        end
        % KL
        [KL, w] = getKL(signal, label,OUTEEG.nMCRSTSfrom);
        
        figure('name','KL-criterion')
        subplot(211)
        plot(krange,w)
        ylabel('Dispersion')
        subplot(212)
        plot(krange(2:end-1),KL)
        ylabel('Krzanovski-Lai criterion')
        
        [~,best_k] = max(KL); % argmax
        
         OUTEEG.nMCRSTS = best_k + OUTEEG.nMCRSTSfrom;
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
            OUTEEG.meanGMD = mean(GMD(signal(1:(size(signal,1)),:),A(:,OUTEEG.idx),size(signal,1)))
        else
            OUTEEG.meanGMD = mean(GMD(signal,A(:,OUTEEG.idx),OUTEEG.nbchan))
        end
        
        if size(signal,2) <= 5000
            show_clusters(signal,OUTEEG.idx,OUTEEG.meanGMD,OUTEEG.nMCRSTS);
        else
            show_clusters(signal(:,1:5000),OUTEEG.idx(1:5000),OUTEEG.meanGMD,OUTEEG.nMCRSTS);
        end
    case 2 % kmeans
        clustering_algorithm = 'kmeans';
        label = zeros(size(signal,2),number_of_ks);
        K = double(fastrbf(signal',rbfsigma)); % double() because otherwise bug in knkmeans in line 22
        %K = doubles(kernelmatrix('poly',signal,signal,0.01,0,4));
        
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
            c
            %show_clusters(signal,label,energy,k);
            nclusters = numel(unique(label(:,k-OUTEEG.nMCRSTSfrom+1)));
            if nclusters < k
                disp(['Couldnt estimate more than ', num2str(nclusters), ' clusters'])
                break
            end
        end
        % KL
        [KL, w] = getKL(signal, label,OUTEEG.nMCRSTSfrom);
        
        figure('name','KL-criterion')
        subplot(211)
        plot(krange,w)
        ylabel('Dispersion')
        subplot(212)
        plot(krange(2:end-1),KL)
        ylabel('Krzanovski-Lai criterion')
        [~,best_k] = max(KL);
        OUTEEG.nMCRSTS = best_k + OUTEEG.nMCRSTSfrom;
        OUTEEG.idx = label(:,best_k+1)';
        A = zeros(OUTEEG.nbchan,OUTEEG.nMCRSTS);
        for iMCRST=1:OUTEEG.nMCRSTS
            A(:,iMCRST) = mean(OUTEEG.data(:,OUTEEG.idx==iMCRST),2);
        end
        OUTEEG.A = A;
        if chronos
            OUTEEG.meanGMD = mean(GMD(signal(1:(size(signal,1)-1),:),A(:,OUTEEG.idx),size(signal,1)-1))
        else
            OUTEEG.meanGMD = mean(GMD(signal,A(:,OUTEEG.idx),OUTEEG.nbchan))
        end
        if size(signal,2) <= 5000
            show_clusters(signal,OUTEEG.idx,0,OUTEEG.nMCRSTS);
        end
        disp('Done with kmeans')
    case 3
        clustering_algorithm = 'agglomerative clustering';
        disp('Running agglomerative clustering...')
        distmat = pdist(signal',@GMD, OUTEEG.nbchan);
        Z = PHA_Clustering(squareform(distmat));
        
        if size(signal,2) < 2000
            figure('name','dendrogram')
            dendrogram(Z);
        end
        
        for k=krange
            label(:,k-OUTEEG.nMCRSTSfrom+1) = cluster(Z,'maxclust',k);
        end
        % KL
        [KL, w] = getKL(signal, label,OUTEEG.nMCRSTSfrom);
        
        figure('name','KL-criterion')
        subplot(211)
        plot(krange,w)
        ylabel('Dispersion')
        subplot(212)
        plot(krange(2:end-1),KL)
        ylabel('Krzanovski-Lai criterion')
        [~,best_k] = max(KL);
        
        OUTEEG.nMCRSTS = best_k + OUTEEG.nMCRSTSfrom;
        OUTEEG.idx = label(:,best_k+1)';
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
        if size(signal,2) <= 5000
            show_clusters(signal,OUTEEG.idx,0,OUTEEG.nMCRSTS);
        end
        disp('Done with agglomerative clustering')
    otherwise
        disp('This is not an option.')
end
disp(['Estimated ', num2str(OUTEEG.nMCRSTS), ' microstates, using ', clustering_algorithm, ' on ', subset, ' of the data']);

% return the string command
% -------------------------
com = sprintf('pop_getmicrostates( %s, %d, [%s] );', inputname(1), int2str(OUTEEG.subset), int2str(OUTEEG.nMCRSTSfrom), int2str(OUTEEG.nMCRSTSto), int2str(OUTEEG.clustering_algorithm));

return;
