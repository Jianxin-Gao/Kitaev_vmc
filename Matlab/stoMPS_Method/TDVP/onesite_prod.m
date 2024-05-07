function T = onesite_prod(T, H_local, H_envi_R, H_envi_L)
T = contract(T, 1, H_envi_L, 2);
T = contract(T, [2 4], H_local, [4 1]);
T = contract(T, [1 3], H_envi_R, [2 3]);
T = permute(T, [1 3 2]);
end