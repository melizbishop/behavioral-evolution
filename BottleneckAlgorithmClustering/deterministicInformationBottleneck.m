function [f,IYT,HT,pY_T,pT] = deterministicInformationBottleneck(pXY,k,f0,beta,tol,maxIter)
    

    if iscell(pXY)
        a = unique(pXY{1}(:));
        b = unique(pXY{2}(:));
        pXY = hist3([pXY{1}(:) pXY{2}(:)],{a,b});
    end
    
    pXY = pXY ./ sum(pXY(:));
    pX = sum(pXY);
    pY_X = bsxfun(@rdivide,pXY,pX);
    
    s = size(pXY);
    N = s(1);
    M = s(2);
    
    if nargin < 3 || isempty(f0)
        f = randi(k,[N 1]);
    end
    
    
    pT = zeros(k,1);
    pY_T = zeros(k,M);
    
    for i=1:k
        pT(i) = sum(pX(f==i));
    end
    idx = pT > 0;
    HT = -sum(pT(idx).*log2(pT(idx)));
    
    for i=1:k
        if pT(i) > 0
            pY_T(i,:) = sum(pXY(f==i,:),1)'./pT(i);
        else
            pY_T(i,:) = 0;
        end
    end
    pYT = bsxfun(@times,pY_T,pT);
    pY = sum(pYT,1);
    temp = pYT.*log2(pYT./(pT*pY));
    IYT = sum(temp(~isnan(temp) & ~isinf(temp)));
    
    endLoop = false;
    n = 1;
    while ~endLoop
        
        previousJ = HT - beta*IYT;
        
        DKLs = findListKLDivergences(pY_X',pY_T);
        fMat = bsxfun(@minus,log2(pT'),beta.*DKLs);
        
        [~,f] = max(fMat,[],2);
        
        for i=1:k
            pT(i) = sum(pX(f==i));
        end
        idx = pT > 0;
        HT = -sum(pT(idx).*log2(pT(idx)));
        
        
        for i=1:k
            if pT(i) > 0
                pY_T(i,:) = sum(pXY(f==i,:),1)'./pT(i);
            else
                pY_T(i,:) = 0;
            end
        end
        pYT = bsxfun(@times,pY_T,pT);
        pY = sum(pYT);
        temp = pYT.*log2(pYT./(pT*pY));
        IYT = sum(temp(~isnan(temp) & ~isinf(temp)));
        
        J = HT - beta*IYT;
        
        if abs(J-previousJ) < tol || n >= maxIter
            break;
        else
            n = n + 1;
        end
        
    end
    
    vals = unique(f);
    if length(vals) < k
        for i=1:length(vals)
            f(f == vals(i)) = i;
        end
        pY_T = pY_T(vals,:);
        pT = pT(vals);
    end
    
    
    
    