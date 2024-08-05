function v = kementrality_old_negative_ldl(G, weights, parallel)
% Computes the "kemeny-based centrality" of the edges of G 
%
% weights (optional) is a vector of weights
%
% Compute the old version of the Kemeny-based centrality used in 
% https://www.hvl.no/globalassets/hvl-internett/arrangement/2022/13sss/522altafini.pdf. 
% Note that this measure can take negative values, and relies on rank-2
% updates when deleting a road disconnects the graph. 

if not(exist('parallel', 'var')) || isempty(parallel)
    parallel = false;
end

n = numnodes(G);
m = numedges(G);

v = nan(m, 1);

coco = conncomp(G);
edges_coco = coco(G.Edges.EndNodes(:,1));
if max(coco) > 1
    %disconnected graph: call this function recursively
    for c = 1:max(coco)
        H = subgraph(G,find(coco==c)); % the 'find' is necessary for this to work on R2017a
        v(edges_coco==c) = kementrality_edges6(H, weights(edges_coco==c));
    end
    return
end

if exist('weights', 'var')
    A = adjacency2(G, weights);
else
    A = adjacency(G);
end

e = ones(n, 1);
w = sparse(n, 1); w(end) = 1;

d = sum(A, 2);

%we cache this transpose since accesses seem slow
Pt = d .\ A; Pt = Pt';

% [L, R] = lu(speye(n) - P + sparse(e)*w');

% blocks of M = D-A+d*w'
M11 = diag(d(1:n-1)) - A(1:end-1, 1:end-1);
M12 = full(-A(1:end-1, end) + d(1:n-1));
M21 = -A(end, 1:end-1);
M22 = full(d(end) - A(end, end) + d(end));
%p = amd(M11); % determined to work better with some experiments
%C = chol(M11(p, p)); % we cannot use this since in some of the larger examples positive
%definiteness seems to be lost numerically
% TODO: unclear if we still need the permutation once we switch to LDL
fprintf('Starting factorization...');
tic
[L, D, p] = ldl(M11, 'vector');
fprintf('Done\n');
toc
% computed solves M11^{-1}(b)
function b = inv_M11(b)
    b(p,:) = L'\(D\(L\(b(p,:))));
end
inv_M11_M12 = inv_M11(M12);
SchurC = M22 - M21*(inv_M11_M12);
fulld = full(d);

% Solves (speye(n) - P + sparse(e)*w')x = b via Cholesky on its
% (1:n-1,1:n-1) submatrix.
% (experimentally, this speeds solves a lot)
function b = solve_via_chol(b)
    b = fulld .* b;
    inv_M11_b1 = inv_M11(b(1:end-1,:));
    b(end,:) = SchurC \ (b(end,:) - M21*inv_M11_b1);
    b(1:end-1,:) = inv_M11_b1 - inv_M11_M12*b(end,:);
end

solve_via_chol_handle = @(b) solve_via_chol(b);

M = G.Edges.EndNodes';

if parallel
    maxWorkers = inf;
else
    maxWorkers = 0;
end
fprintf('Started centrality computation...\n')
tic
checkpoints = unique(floor(linspace(1,m,50)));
fprintf('Each character # is printed when approximately 1/%d of the computation is completed:\n',length(checkpoints));
parfor (k = 1:size(M,2), maxWorkers)
%for k = 1:size(M,2)  % one may want to switch to a non-parallel for for profiling
    if ismember(k, checkpoints)
        fprintf('#\n');
    end

    ij = M(:,k);
    
    %
    % we set up U, V so that P + U*V' is the random walk matrix
    % of the graph with edge ij removed.
    %
    
    U = zeros(n, 2);
    U(ij(1), 1) = 1;
    U(ij(2), 2) = 1;
    Pij = full(Pt(:,ij)');
    newrows = Pij;
    newrows(1, ij(2)) = 0;
    newrows(2, ij(1)) = 0;
    diagsum = sum(newrows,2);
    for h = 1:2 % add a self-edge in case a row becomes empty after the removal
        if diagsum(h) == 0
            diagsum(h) = 1;
            newrows(h, ij(h)) = 1;
        end
    end
    newrows = sum(newrows,2) .\ newrows;
    V = newrows - Pij;

    %
    % Sherman-Morrison-based update
    %
    
%    invU = R  \ (L \ U);
    invU = solve_via_chol_handle(U);
    % solving systems with column vectors seems faster than systems with
    % row vectors, so we compute VA^{-1}A^{-1}U via A^{-1}A^{-1}U
    % invinvU = R \ (L \ invU);
    invinvU = solve_via_chol_handle(invU);
    core = eye(2) - V * invU;
    H = rmedge(G, k);
    c = conncomp(H);
    if max(c) == 1
        v(k) = trace(V * invinvU / core);
    else % our update disconnected the graph; we need to make another update to fix this
        % the update is in the form kervec*z';
        kervec = double(c == 1)';
        z = [ones(n-2,1); 1; 0];
        z = z / (z'*kervec); %normalizes the update vector
        %invU = [invU,  R \ (L \ kervec)];
        invU = [invU, solve_via_chol_handle(kervec)];
        %invinvU = [invinvU R \ (L \ invU(:, end))];
        invinvU = [invinvU solve_via_chol_handle(invU(:, end))];
        V = [V; z'];
        core = eye(3) - V * invU;
        v(k) = trace(V * invinvU / core) + 1;
    end
end
timeparallel = toc;
fprintf('Done! Cpu time %f s.\n', timeparallel)

end