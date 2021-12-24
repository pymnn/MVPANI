function ranking = llcfs(X)
% Feature Selection and Kernel Learning for Local Learning-Based Clustering
% Input
%	X: nSmp * nDim
%	param, a struct of parameters
%		nClusters, the number of clusters
%		k, the size of knn
%		beta, the regularization parameter
% Output
%	Y: nSmp * nClusters
%	tao: nDim * 1
%	
% [1] Feature Selection and Kernel Learning for Local Learning-Based Clustering, PAMI-2011
%  Version 5.0 August 2017
%  Support: Giorgio Roffo
%  E-mail: giorgio.roffo@glasgow.ac.uk

%  If you use our toolbox please cite our supporting papers:
% 
%  BibTex
%  ------------------------------------------------------------------------
% @InProceedings{RoffoICCV17, 
% author={Giorgio Roffo and Simone Melzi and Umberto Castellani and Alessandro Vinciarelli}, 
% booktitle={2017 IEEE International Conference on Computer Vision (ICCV)}, 
% title={Infinite Latent Feature Selection: A Probabilistic Latent Graph-Based Ranking Approach}, 
% year={2017}, 
% month={Oct}}
%  ------------------------------------------------------------------------
% @Inbook{Roffo2017,
% author="Roffo, Giorgio
% and Melzi, Simone",
% editor="Appice, Annalisa
% and Ceci, Michelangelo
% and Loglisci, Corrado
% and Masciari, Elio
% and Ra{\'{s}}, Zbigniew W.",
% title="Ranking to Learn:",
% bookTitle="New Frontiers in Mining Complex Patterns: 5th International Workshop, NFMCP 2016, Held in Conjunction with ECML-PKDD 2016, Riva del Garda, Italy, September 19, 2016, Revised Selected Papers",
% year="2017",
% publisher="Springer International Publishing",
% address="Cham",
% pages="19--35",
% isbn="978-3-319-61461-8",
% doi="10.1007/978-3-319-61461-8_2",
% url="https://doi.org/10.1007/978-3-319-61461-8_2"
% }
%  ------------------------------------------------------------------------
% @InProceedings{RoffoICCV15, 
% author={G. Roffo and S. Melzi and M. Cristani}, 
% booktitle={2015 IEEE International Conference on Computer Vision (ICCV)}, 
% title={Infinite Feature Selection}, 
% year={2015}, 
% pages={4202-4210}, 
% doi={10.1109/ICCV.2015.478}, 
% month={Dec}}
%  ------------------------------------------------------------------------
%===================setup=======================
knnCandi = 5;
graphTypeCandi = [2];
betaCandidates = 10.^-1;
paramcell = llcfs_build_param(2, knnCandi, betaCandidates, graphTypeCandi );
%===============================================
param = paramcell{1};
X = double(X);

if isfield(param, 'nClusters')
    c = param.nClusters;
end

k = 30;
if isfield(param, 'k')
    k = param.k;
end

beta = 1;
if isfield(param, 'beta')
    beta = param.beta;
end

kType = 1;
if isfield(param, 'kType')
    kType = param.kType;
end

maxiter = 2;
if isfield(param, 'maxiter')
    maxiter = param.maxiter;
end

epsilon = 1e-5;
if isfield(param, 'epsilon')
    epsilon = param.epsilon;
end

isTao = 0;
epsilon_tao = 1e-5;
[n, d] = size(X);


% convergence by maxiter
isMaxiter = 1;
if maxiter > 0
    isMaxiter = 1;
end

% convergence by epsilon
isEpsilon = 0;
if isEpsilon > 0
    isEpsilon = 1;
end

tao = ones(d,1) / d;

objHistory = [];
iter = 0;
while true
    
    wX = bsxfun(@times, X, sqrt(max(tao, eps))' );
    wX2 = bsxfun(@times, X, max(tao, eps)' );
    wK = wX * wX';
    % k-mutual neighbors re-computation using weighted features
    switch kType
        case 1
            W = SimGraph_NearestNeighbors(wX', k, 2, 0);
            [idx, jdx, ~] = find(W);
            kIdx = cell(n, 1);
            nz = length(idx);
            for ii = 1:nz
            	kIdx{jdx(ii)} = [kIdx{jdx(ii)}, idx(ii)];
            end
        case 2
            if isempty(which('knnsearch'))
                disp('The funcion knnsearch in stat toolbox is not found');
            else
                [kIdx, ~] = knnsearch(wX, wX, 'k', min(n, k + 1) );
                kIdx = kIdx(:, 2:end);
                kIdx = mat2cell(kIdx, ones(n, 1), size(kIdx, 2));
            end
        otherwise
            disp('');
    end
    
    % construct A for laplacian
    A = zeros(n);
    wA = cell(n,1);% pre storage for w computation
    for i = 1:n
        lidx = kIdx{i};
        ni = length(lidx);
        if ni > 1
            Ki = wK(lidx, lidx);
            ki = wK(i, lidx);
            Hi = eye(ni) - ones(ni, ni) / ni;
            Ii = eye(ni);
            Iib = Ii / beta;
            Ai = Hi * Ki * Hi;
            Ai = (Ai + Iib) \ Ai;
            Ai = Hi - Hi * Ai;
            Ai = Ai * beta;
            wA{i} = wX2(lidx, :)' * Ai; % EQ 15
            Ai = (ki - sum(Ki) / ni) * Ai;
            Ai = Ai + ones(1, ni) / ni;
            A(i, lidx) = Ai;
        end
    end
    
    % construct laplacian for local learning
    M = eye(n) - A;
    M = M' * M;
    M(isnan(M)) = 0;
	M(isinf(M)) = 0;
	
    % first c eigenvectors corresponding to the first c smallest eigenvalues
    M = (M + M') / 2;
    [Y, eigval] = eig(M);
    eigval = diag(eigval);
    [eigval, eigidx] = sort(eigval, 'ascend');
	eigval = eigval(eigidx(1:c));
    Y = Y(:, eigidx(1:c));
    
    objHistory = [objHistory; sum(eigval)];%#ok
    
	
    % compute wc to compute tao
    tao_old = tao;
	
	tao = zeros(d, 1);
    for i = 1:n
        lidx = kIdx{i};
        ni = length(lidx);
        if ni > 1
            wi = wA{i} * Y(lidx,:);
            tao = sum(wi.^2, 2) + tao;
        end
    end
	tao = sqrt(tao);
    tao = tao / sum(tao);
    
    % check the convergence
    iter = iter + 1;
    if isEpsilon && iter > 1
        if abs(objHistory(end-1) - objHistory(end)) < epsilon
            break;
        end
    end
	if isTao && sum(abs(tao_old - tao)) < epsilon_tao
		break;
	end
    if isMaxiter && iter == maxiter
        break;
    end
end

[~, ranking] = sort(tao, 'descend');



end

