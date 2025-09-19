function label = DGGCE(BPi, nSub, nSubBase, seed, nCluster)
cntTcutKmReps = 10; 
maxTcutKmIters = 100; 
Bs = cell(1, nSub); 

%*********************************
% 生成随机标签组
%*********************************
[~, BPi_subspaces] = generate_random_sublabelgroups(BPi, nSub, nSubBase, seed);


for i = 1:nSub
    Hc = compute_Hc(BPi_subspaces{i});
    Bs{i} = GB_USC_v1(Hc);
end



HcB = cell2mat(Bs);


Order = 5;
P = cell(1, Order);
options = [];
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'HeatKernel';

Hc_tmp = zscore(HcB);
Ds = EuDist2(Hc_tmp', Hc_tmp');
t = max(median(Ds, 'all'), eps);

options.t = t;
options.bSelfConnected = 1;
nSmp_S = size(HcB, 2);
S = constructW(HcB', options);
D = diag(sum(S, 2));
P_1 = 0.5 * (eye(nSmp_S) + D^(-0.5) * S * D^(-0.5));


P{1} = P_1;
for k_idx = 2:Order
    P{k_idx} = P{k_idx - 1} * P_1;
end
W = ones(1, Order) / Order;
G = zeros(nSmp_S);
for k_idx = 1:Order
    P_norm = norm(P{k_idx}, 'fro');
    if P_norm > 0
        P{k_idx} = P{k_idx} / P_norm;
    end
    G = G + W(k_idx).*P{k_idx};
end
HcB_new  = HcB * G;

[label,~] = Tcut_for_bipartite_graph(HcB_new,nCluster,maxTcutKmIters,cntTcutKmReps);


end

function [labels,evec] = Tcut_for_bipartite_graph(B,Nseg,maxKmIters,cntReps)
if nargin < 4
    cntReps = 3;
end
if nargin < 3
    maxKmIters = 100;
end

[Nx,Ny] = size(B);
if Ny < Nseg
    error('Need more columns!');
end

dx = sum(B,2);
dx(dx==0) = 1e-10;
Dx = sparse(1:Nx,1:Nx,1./dx); clear dx
Wy = B'*Dx*B;

d = sum(Wy,2);
D = sparse(1:Ny,1:Ny,1./sqrt(d)); clear d
nWy = D*Wy*D; clear Wy
nWy = (nWy+nWy')/2;


[evec,eval] = eig(full(nWy)); clear nWy
[~,idx] = sort(diag(eval),'descend');
Ncut_evec = D*evec(:,idx(1:Nseg)); clear D

evec = Dx * B * Ncut_evec; clear B Dx Ncut_evec

evec = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );

% k-means
labels = kmeans(evec,Nseg,'MaxIter',maxKmIters,'Replicates',cntReps);
end