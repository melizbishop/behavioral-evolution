function [clusterings,HTs,IYTs,numClusters,betas,clusterChoices] = ...
            run_DIB(X,Y,N,minClusters,maxClusters,minLogBeta,maxLogBeta)

    %Inputs:
    %X -> L x 1 array of first variable
    %Y -> L x 1 array of second variable
    %N -> Number of iterations (default = 10000)
    %minClusters -> minimum number of clusters (default = 2)
    %maxClusters -> maximum number of clusters (default = 30)
    %minLogBeta -> minimum value of log10(beta) (default = -1)
    %maxLogBeta -> maximum value of log10(beta) (default = 4)
    %
    %
    %Outputs:
    %clusterings -> R x 1 cell array of pareto front clusterings
    %HTs -> R x 1 array of entropies for each clustering
    %IYTs -> R x 1 array of I(Y;T) for each clustering
    %numClusters -> R x 1 array of number of clusters for each clustering
    %betas -> R x 1 array of betas for each clustering
    %clusterChoices -> R x 1 binary array. True if clustering 
    %               has the largest I(Y;T) for a given number of clusters 
    
    
    if nargin < 3 || isempty(N)
        N = 10000;
    end
    
    if nargin < 4 || isempty(minClusters)
        minClusters = 2;
    end
    
    if nargin < 5 || isempty(maxClusters)
        maxClusters = 30;
    end
    
    if nargin < 6 || isempty(minLogBeta)
        minLogBeta = -1;
    end
    
    if nargin < 7 || isempty(maxLogBeta)
        maxLogBeta = 4;
    end
    
    readout = 100;

    betas = zeros(N,1);
    numClusters = zeros(N,1);
    clusterings = cell(N,1);
    IYTs = zeros(N,1);
    HTs = zeros(N,1);
    
    a = unique(X(:));
    b = unique(Y(:));
    pXY = hist3([X(:) Y(:)],{a,b});
    
    parfor i=1:N
        betas(i) = 10^(minLogBeta + (maxLogBeta-minLogBeta)*rand());
        k = minClusters + randi(maxClusters - minClusters) - 1;
        [clusterings{i},IYTs(i),HTs(i),~,~] = ...
            deterministicInformationBottleneck(pXY,k,[],betas(i),1e-6,1000);
        numClusters(i) = length(unique(clusterings{i}));
        if mod(i,readout) == 0
            fprintf(1,'Calculating for Iteration #%6i out of %6i\n',i,N);
        end
    end
    
    idx = findParetoFront([-HTs IYTs]);
    clusterings = clusterings(idx);
    IYTs = IYTs(idx);
    HTs = HTs(idx);
    numClusters = numClusters(idx);
    
    [~,sortIdx] = sort(IYTs);
    clusterings = clusterings(sortIdx);
    IYTs = IYTs(sortIdx);
    HTs = HTs(sortIdx);
    numClusters = numClusters(sortIdx);
    
    idx = [1;find(diff(IYTs) > 1e-10)+1];
    clusterings = clusterings(idx);
    IYTs = IYTs(idx);
    HTs = HTs(idx);
    numClusters = numClusters(idx);
    
    
    clusterValues = unique(numClusters);
    clusterChoices = false(size(numClusters));
    for i=1:length(clusterValues)
        idx = find(numClusters == clusterValues(i),1,'last');
        clusterChoices(idx) = true;
    end
    
    
    