function [kemder,ij] = kderivative(g, weights, parallel, relative)
% function kemder = kderivative(g, weights, parallel, relative)
% Computes the centralities of all the edges of a graph 
% based on the derivative of the Kemeny constant
% In this version, the permutation  made by the function 'amd' is
% not applied inside the loop, and the ordering of the indices
% ij is independent of the graph / adjacency matrix input
% In Input:
%    -- g: either the graph or the adjacency matrix of the graph
%    -- weights: the vector of weights of the graph. If missing
%       or empty [], it is assumed that all the weights are 1
%    -- parallel: if true, "parfor" is used instead of "for"
%       (default true)
%    -- relative: if true the relative kementrality is computed
%       (default true)
% In Output:
%    --  kemder: the vector of the centralities, i.e., kemder(l)
%       is the centrality of the edge {i(l),j(l)} where 
%       i(l)>j(l) and i(l)<i(l+1), j(l)<j(l+1).
%    -- ij the mx3 matrix with columns i,j,kem such that
%       kem(l) is the kementrality of the edge {i(l),j(l)}
%       where i(l)>j(l) and i(l)<=i(l+1), j(l)<=j(l+1).
%       and m is the number of edges.
  
%%% Input analysis
   debug = true;
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
   if not (exist('parallel','var'))
      parallel = true;
   end
   if not (exist('weights','var'))
      weights = [];
   end
   if not (exist('relative','var'))
      relative = true;
   end   
   reg=eps;
   
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
   kemder = zeros(m,1);
  
%  Preprocessing
   z = spdiags(d*(1+reg),0,n,n)-a;  
   z(n,n)=z(n,n)+1;
   L = chol(z)';    %%%% 1
   dt = L\d;        %%%% 2
   al = norm(d,1); sal = sqrt(al); ph = sal/dt(n); 
   ps = -(al+norm(dt)^2)/dt(n)^2;
   pretime = toc;
   fprintf('preprocessing time=%d\n',pretime)
   if debug
       dL=full(diag(L)); 
       mdL=min(dL); 
       EdL=dL(end); 
       MdL=max(dL);
       fprintf('Diag entries of L:  [min,last,max]=[%d,%d,%d]\n',mdL,EdL,MdL);
   end
   
% Main computation
   fprintf('Each character # is printed when approximately 1/50 of the computation is completed:\n');
   checkpoints = unique(floor(linspace(1,m,50)));
   tic
   parfor (l=1:m, workers)
   %   for k = 1:m  % one may want to switch to a non-parallel for for profiling
       if ismember(l, checkpoints)
         fprintf('#\n');
       end
       p = i(l); q = j(l); 
       w=zeros(n,1); w(q) = 1; w(p) = -1;
       f = L\w;                                                        
       ga = ph*f(n); de = ph*(dt'*f)/sal+ps*f(n);
       den = zeros(n,1);den(n)=de;
       g = f-(ga/sal)*dt-den;                                 
       y = L'\g;                                                      
       kemder(l)=sum(d.*y.^2);                          
       if ~relative
           kemder(l)=kemder(l)*weights(l);
       end
   end
%  permute output consistently with the original data   
   a = sparse(i,j,kemder,n,n) ;a = a+a';
   a(per,per) = a;
   [ia,ja,wa] = find(a);
   i = ia(ia>ja); j = ja(ia>ja); kemder = wa(ia>ja);
   ij = [j,i,kemder];
   comptime = toc;
   fprintf('computation time = %d\n',comptime);
end

