function [Rslt_MPS] = Canonicalization(MPS)
% tic
cell_size = size(MPS);
n = cell_size(2);
Rslt_MPS = cell(1, n);

[U, S, V] = svd(MPS{n}, 'econ');
Rslt_MPS{n} = V';
% fprintf('Canonicalization Step: %d/%d', 1, n)
% toc

for i = (n-1):-1:2
%     tic
    Q_size = size(Rslt_MPS{i+1});
    MPS{i} = contract(MPS{i}, 2, U*S, 1);
    MPS{i} = permute(MPS{i}, [1 3 2]);
    T_size = size(MPS{i});
    T_matrix = reshape(MPS{i}, T_size(1), T_size(2)*T_size(3));
    [U, S, V] = svd(T_matrix, 'econ');
    Rslt_MPS{i} = reshape(V', [], Q_size(1), T_size(3));
%     fprintf('Canonicalization Step: %d/%d', (n+1-i), n)
%     toc
end

% tic
MPS{1} = contract(MPS{1}, 1, U*S, 1);
Rslt_MPS{1} = permute(MPS{1}, [2 1]);
% fprintf('Canonicalization Step: %d/%d', n, n)
% toc
end