function T = twosite_prod(T, direct, H_local1, H_local2, H_envi_R, H_envi_L)
switch direct
    case 'fL2R'
        T = contract(T, 1, H_local1, 3);
        T = contract(T, [2 3], H_local2, [4 1]);
        T = contract(T, [1 3], H_envi_R, [2 3]);
        T = permute(T, [1 3 2]);
    case 'L2R'
        T = contract(T, 1, H_envi_L, 2);
        T = contract(T, [1 5], H_local1, [4 1]);
        T = contract(T, [2 4], H_local2, [4 1]);
        T = contract(T, [1 4], H_envi_R, [2 3]);
        T = permute(T, [1 2 4 3]);
    case 'lL2R'
        T = contract(T, 3, H_local2, 3);
        T = contract(T, [2 3], H_local1, [4 2]);
        T = contract(T, [1 3], H_envi_L, [2 3]);
        T = permute(T, [3 2 1]);    
end