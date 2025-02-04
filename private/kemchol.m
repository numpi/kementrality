function [kem, ij] = kemchol(g, reg, weights, parallel)
% function kem = kemchol(g, reg, weights, parallel)
% Computes the Kemeny-based centrality of all the edges of a
% graph 
% In this version, the permutation  made by the function 'amd' is
% not applied inside the loop, and the ordering of the indices
% ij is independent of the the graph / adjacency matrix input
% In Input:
%    -- g: either the graph or the adjacency matrix of the graph
%    -- reg: regularization parameter (default reg = 1.e-6)
%    -- weights: the vector of weights of the graph. If missing
%       or empty [], it is assumed that all the weights are 1
%    -- parallel: if true, "parfor" is used instead of "for"
%       (default true)
% In Output:
%    -- kem: the vector of the centralities, i.e., kem(l)
%       is the centrality of the edge {i(l),j(l)} where 
%       i(l)>j(l) and i(l)<=i(l+1), j(l)<=j(l+1).
%   -- ij the mx3 matrix with columns i,j,kem such that
%       kem(l) is the kementrality of the edge {i(l),j(l)}
%       where i(l)>j(l) and i(l)<=i(l+1), j(l)<=j(l+1).
%       and m is the number of edges.
  
%%% Input analysis
   tic
   if ~isnumeric(g)
      % graph as first input (optionally with weights)
      if not(exist('weights', 'var')) || isempty(weights)
        a = adjacency(g);
      else
        a = adjacency(g, weights);
      end
      clear g;
   else
      % adjacency matrix as first input
      a = g;
      % check weights
      if exist('weights', 'var') && not(isempty(weights))
        error('Specifying weights is not supported with an adjacency matrix as first argument, only with a graph')
      end
      clear g;
   end
   if not (exist('reg','var'))
      reg = 1.e-6;
   end
   if not (exist('parallel','var'))
      parallel = true;
   end
   if not (exist('weights','var'))
     weights = [];
   end
   
%%% Parallel setting   
   if parallel
        workers = inf;
   else
        workers=0;
   end
   
%%% Permute rows and columns of a
   per=amd(a);
   a=a(per,per);

%%% Extract and check information
   [ia, ja, wa] = find(a);
   if any(ia==ja)
      error('the given adjacency matrix has self-loops; this is not supported for now');
   end
   i = ia(ia>ja); j = ja(ia>ja); weights = wa(ia>ja);
   n = size(a, 1); m = length(weights);
   clear ia ja wa;

%%% Start computing
   d = sum(a,2);
   kem = zeros(m,1);
   sd = sum(d);
   ij=[i,j]'; 
   
%  Preprocessing
   T = spdiags(d*(1+reg),0,n,n)-a;  
   R = chol(T);
   td = R' \ d;
   td = R \ td;
   dtd = sum(d .* td);
   sdtd = dtd + sd;
   fprintf('Started centrality computation...\n')
   tic
   checkpoints = unique(floor(linspace(1,m,50)));
   fprintf('Each character # is printed when approximately 1/%d of the computation is completed:\n',length(checkpoints));
   parfor (k = 1:m, workers)
%   for k = 1:m  % one may want to switch to a non-parallel for for profiling
     if ismember(k, checkpoints)
        fprintf('#\n');
     end
     ij_slice = ij(:,k);
     ii = ij_slice(1);
     jj = ij_slice(2);
     Aij = weights(k);
     v = sparse([ii jj], 1, [1 -1], n, 1);
     v=full(v);
     w = R' \ v;
     w = R \ w;
     dw = sum(d .* w);
     x = w - (dw/sdtd) * td;
     alpha = Aij * (x(ii)-x(jj));
     beta = Aij * sum(x .* x .* d);
     ce = beta / (1-alpha);
     if ce > 0.5/reg
       ce = abs(ce - 1/reg);
     end
     kem(k) = ce;
   end
   timeparallel = toc;
   fprintf('Done! Cpu time %f s.\n', timeparallel)

%  permute output consistently with the original data  
   a = sparse(i,j,kem,n,n);
   a = a+a';
   a(per,per) = a;
   [ia,ja,wa] = find(a);
   i = ia(ia>ja); j = ja(ia>ja); kem = wa(ia>ja);
   ij=[j,i,kem];
   comptime = toc;
   fprintf('computation time = %d\n',comptime);
end
   



