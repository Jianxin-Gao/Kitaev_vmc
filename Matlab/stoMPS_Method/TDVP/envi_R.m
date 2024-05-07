function [H_envi] = envi_R(MPS_cell, H)
n = length(MPS_cell);
H_envi = cell(1, n);
H_envi{n} = contract(H{n}, 3, MPS_cell{n}, 2);
H_envi{n} = contract(H_envi{n}, 2, conj(MPS_cell{n}), 2);
H_envi{n} = permute(H_envi{n}, [3 2 1]);
for i = n-1:-1:3
    H_envi{i} = contract(H_envi{i+1}, 2, MPS_cell{i}, 2);
    H_envi{i} = contract(H_envi{i}, [2 4], H{i}, [2 4]);
    H_envi{i} = contract(H_envi{i}, [1 4], conj(MPS_cell{i}), [2 3]);
    H_envi{i} = permute(H_envi{i}, [3 1 2]);
end
end