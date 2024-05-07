function [H_envi_R, MPS_cell_rev, H_rev] = reverse(H_envi_L, MPS_cell, H)
n = length(MPS_cell);
H_envi_R = flip(H_envi_L);
MPS_cell_rev = flip(MPS_cell);
H_rev = flip(H);
for i = 2:n-1
    MPS_cell_rev{i} = permute(MPS_cell_rev{i}, [2 1 3]);
    H_rev{i} = permute(H_rev{i}, [2 1 3 4]);
end
end