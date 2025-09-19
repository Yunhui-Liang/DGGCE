function [Ks] = X2Ks_large_20250724_same(X1, X2, k, blockSize)

% 输入：
% X1, X2：两个数据集，X1 (n1 x d), X2 (n2 x d)
% k: 选取的 top-k 样本数
% blockSize: 用于分块计算的大小

[nSmp1, d] = size(X1);
[nSmp2, ~] = size(X2);

nKernel = 12;
Ks = cell(1, nKernel);

if ~exist('blockSize', 'var')
    blockSize = 1000;
end
blockSize = min(blockSize, nSmp1); 
blockSize = min(blockSize, nSmp2);
numBlocks = ceil(nSmp1 / blockSize);
rowIdx_s_cell = cell(numBlocks, 1);
colIdx_s_cell = cell(numBlocks, 1);
val_s_cell = cell(numBlocks, 1);
rowIdx_d_cell = cell(numBlocks, 1);
colIdx_d_cell = cell(numBlocks, 1);
val_d_cell = cell(numBlocks, 1);

t_k = 0;
x_norms1 = sum(X1.^2, 2); 
x_norms2 = sum(X2.^2, 2);
x_min = min([x_norms1; x_norms2]);
x_max = max([x_norms1; x_norms2]);

s = zeros(1, numBlocks); 

for iBlock = 1:numBlocks 
    blockStart = (iBlock - 1) * blockSize + 1; 
    blockEnd = min(iBlock * blockSize, nSmp1);
    blockSizeCurr = blockEnd - blockStart + 1; 
    
    t1 = tic;
    % 计算X1和X2之间的核值
    s_block = X1(blockStart:blockEnd, :) * X2';
    
    if x_min > 0.9999 && x_max < 1.001   
        s(iBlock) = 2 * nSmp1 * blockSizeCurr - 2 * sum(sum(s_block));
    else  
        d_block = bsxfun(@plus, -2 * s_block, x_norms1(blockStart:blockEnd));
        d_block = bsxfun(@plus, d_block, x_norms2'); 
        s(iBlock) = sum(sum(d_block));
    end
    t_k = t_k + toc(t1); 
    if nSmp1 == nSmp2
        % Subtract diagonal values for block
        s_block(:, blockStart:blockEnd) = s_block(:, blockStart:blockEnd) - 10^8 * eye(blockSizeCurr);
    end

    % Sort and select top-k
    [K_block_value, Idx] = sort(s_block, 2, 'descend'); 
    if length(Idx(1,:)) <= k
        k = length(Idx(1,:));
    end
    Idx = Idx(:, 1:k);  
    K_block_value = K_block_value(:, 1:k); 

    rowIdx2 = repmat((1:blockSizeCurr)', k, 1); 
    rowIdx_s_cell{iBlock} = rowIdx2(:) + blockStart - 1;
    colIdx_s_cell{iBlock} = Idx(:);
    val_s_cell{iBlock} = K_block_value(:); 

    % Handle diagonal subtraction and sorting
    if nSmp1 == nSmp2
        d_block(:, blockStart:blockEnd) = d_block(:, blockStart:blockEnd) + 10^8 * eye(blockSizeCurr);
    end
    [D_block_value, Idx] = sort(d_block, 2, 'ascend');
    Idx = Idx(:, 1:k); 
    D_block_value = D_block_value(:, 1:k);

    rowIdx2 = repmat((1:blockSizeCurr)', k, 1);
    rowIdx_d_cell{iBlock} = rowIdx2(:) + blockStart - 1;
    colIdx_d_cell{iBlock} = Idx(:); 
    val_d_cell{iBlock} = D_block_value(:); 
end

% Combine indices and values for K
rowIdx = cell2mat(rowIdx_s_cell);
colIdx = cell2mat(colIdx_s_cell);
val_s = cell2mat(val_s_cell);
Ks{1} = create_K(rowIdx, colIdx, val_s, nSmp1, nSmp2);
val = (1 + val_s).^2;
Ks{2} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);
val = (1 + val_s).^4;
Ks{3} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);
val = val_s.^2;
Ks{4} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);
val = val_s.^4;
Ks{5} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

% Kernel for D (distance)
rowIdx = cell2mat(rowIdx_d_cell);
colIdx = cell2mat(colIdx_d_cell);
val_d = cell2mat(val_d_cell);

s2 = sum(s) / nSmp1^2;
val = exp(-val_d / (2^-3 * s2));
Ks{6} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

% Gaussian kernels with different bandwidths
val = exp(-val_d / (2^-2 * s2));
Ks{7} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

val = exp(-val_d / (2^-1 * s2));
Ks{8} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

val = exp(-val_d / (2^0 * s2));
Ks{9} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

val = exp(-val_d / (2^1 * s2));
Ks{10} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

val = exp(-val_d / (2^2 * s2));
Ks{11} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

val = exp(-val_d / (2^3 * s2));
Ks{12} = create_K(rowIdx, colIdx, val, nSmp1, nSmp2);

end

% Helper function to create the kernel
function K = create_K(rowIdx, colIdx, val, nSmp1, nSmp2)
K = sparse(rowIdx, colIdx, val, nSmp1, nSmp2);  
% K = full(K);
% K = bsxfun(@rdivide, K, sum(K, 2) + eps);
if nSmp1 == nSmp2
    K = K + speye(nSmp1);
    K = 0.5 * K + 0.5 * K';
end
% DSsym = 1 ./ sqrt(sum(K, 1) + eps);
% K = bsxfun(@times, K, DSsym);
% K = bsxfun(@times, K, DSsym');
% K = K + speye(nSmp);
% K = 0.5 * K;
% K = 0.5 * K + 0.5 * K';
end 