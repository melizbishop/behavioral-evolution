function [P,dists,infos,Fs,endDists,endInfos] = ...
                    information_cluster_pairwise_distances(D,k,T,maxIter,epsilon,initialP,displayOff)

   %Inputs:
   %D -> NxN distance matrix
   %k -> # of clusters
   %T -> temperature (warning: all sorts of bugs happen if you make T too
   %        small) (default = 1)
   %maxIter -> Maximum # of iterations (default = 1000)
   %epsilon -> relative convergence critereon (default = 1e-5)
   %initialP -> initial condition for P ([] if none, default = [])
   %displayOff -> toggles off output during run (default = false)
   
   %Outputs:
   %P -> Nxk probablility matrix of each data point (row) belonging to a
   %        cluster (column)
   %dists -> average within cluster distance as a function of iteration
   %infos -> average information as a function of iteration
   %Fs -> Cost function value as a function of iteration
   %endDists -> within cluster distances at the end of the calculation
   %endInfos -> within cluster contributions to the information at the 
   %        end of the calculation
   
   
   
   if nargin < 3 || isempty(T)
       T = 1;
   end
   
   if nargin < 4 || isempty(maxIter)
       maxIter = 1000;
   end
   
   if nargin < 4 || isempty(epsilon)
       epsilon = 1e-5;
   end
   
   if nargin < 7 || isempty(displayOff)
       displayOff = false;
   end
   
   MINVALUE = 1e-50;
   
   N = length(D(:,1));
   %D2 = .5*(D + D');
   diffNorms = zeros(maxIter,1);
   maxDiffNorms = zeros(maxIter,1);
   dists = zeros(maxIter,1);
   infos = zeros(maxIter,1);
   Fs = zeros(maxIter,1);
   
   if nargin < 6 || isempty(initialP)
       
       P = rand(N,k);
       P = P ./ repmat(sum(P,2),1,k);
       
   else
       P = initialP;
       P = P ./ repmat(sum(P,2),1,k);
   end
   
   simVals = zeros(N,k);
   Zs = zeros(N,k);
   currentInfos = zeros(1,k);
   
   count = 1;
   test = true;
   
   
   while test
       
       if ~displayOff
           if count > 1
               fprintf(1,'Iteration #%5i\tDists = %8f\tInfos = %8f\tF =%8f\n',count,dists(count-1),infos(count-1),Fs(count-1));
           else
               fprintf(1,'Iteration #%5i\n',count);
           end
       end
       
       Z = sum(P);
       Z(Z < MINVALUE) = MINVALUE;
       
       tempSimVals = diag(P'*D*P) ./ Z'.^2;
       
       for i=1:k
           q = P(:,i);
           if sum(q) > 0
               z = q(q > 0);
               currentInfos(i) = sum(z.*(log(z) - log(Z(i)) + log(N))./log(2))/N;
           else
               currentInfos(i) = 0;
           end
           simVals(:,i) = tempSimVals(i);
       end
       
       
       for i=1:k
           Zs(:,i) = Z(i);
       end
       
       simMinusVals = (D*P) ./ Zs;
       PC = repmat(mean(P),N,1);
       
       
       dists(count) = sum(PC(1,:).*simVals(1,:));
       infos(count) = sum(currentInfos);
       
       Pnew = PC.*exp((simVals - 2*simMinusVals)./T);
       sum_Pnew = sum(Pnew,2);
       sum_Pnew(sum_Pnew < MINVALUE) = MINVALUE;
       Pnew = bsxfun(@rdivide,Pnew,sum_Pnew);
       
       
       q = Pnew - P;
       diffNorms(count) = norm(q,1);
       maxDiffNorms(count) = max(abs(q(:)));
       Fs(count) = dists(count) + T*infos(count);
       
       
       
       P = Pnew;
       
       if count == 1
           
           count = count + 1;
           
       else
           
           if (Fs(count-1) - Fs(count))/Fs(count-1) < epsilon || count == maxIter
               test = false;
           else
               count = count + 1;
           end
           
       end
       
   end
   
   dists = dists(1:count);
   infos = infos(1:count);
   Fs = Fs(1:count);
   
   endDists = simVals(1,:);
   endInfos = currentInfos;
   
   