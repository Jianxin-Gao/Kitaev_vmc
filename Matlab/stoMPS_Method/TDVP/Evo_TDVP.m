function [nMPS1, nMPS2, H_envi, log_S_norm, t_error, Se, Kry_K] = Evo_TDVP(MPS1, MPS2, tau, D, direct, H_envi_R, H_envi_L, H_local1, H_local2)
switch direct
    case 'fL2R'
        MPS_temp = contract(MPS1, 1, MPS2, 1);
        %         MPS_temp = contract(MPS_temp, 1, H_local1, 3);
        %         MPS_temp = contract(MPS_temp, [2 3], H_local2, [4 1]);
        %         MPS_temp = contract(MPS_temp, [1 3], H_envi_R, [2 3]);
        %         MPS_temp = permute(MPS_temp, [1 3 2]);
        
        fun_prod = @(T) twosite_prod(T, 'fL2R', H_local1, H_local2, H_envi_R, H_envi_L);
        fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
        fun_norm = @(T) tensor_norm(T);
        [MPS_temp, Info] = LanczosExp(fun_prod, MPS_temp, -tau, 'innerprod', fun_inner, 'norm', fun_norm);
        Kry_K1 = Info.K;
        Tsize = size(MPS_temp);
        MPS_mat = reshape(MPS_temp, Tsize(1), Tsize(2)*Tsize(3));
        [U, S, V, t_error, Se] = svdT(MPS_mat, 'Nkeep', D, 'epsilon', 1e-12);
        log_S_norm = log(norm(S));
        S = S/norm(S);
        S_size = size(S);
        nMPS1 = permute(U, [2 1]);
        nMPS2 = reshape(diag(S)*V, S_size(1), [], 2);
        
        H_envi = contract(H_local1, 3, nMPS1, 2);
        H_envi = contract(H_envi, 2, conj(nMPS1), 2);
        H_envi = permute(H_envi, [3 2 1]);
        nMPS_temp = nMPS2;
        %         nMPS_temp = contract(nMPS2, 1, H_envi, 2);
        %         nMPS_temp = contract(nMPS_temp, [2 4], H_local2, [4 1]);
        %         nMPS_temp = contract(nMPS_temp, [1 3], H_envi_R, [2 3]);
        %         nMPS_temp = permute(nMPS_temp, [1 3 2]);
        
        fun_prod = @(T) onesite_prod(T, H_local2, H_envi_R, H_envi);
        fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
        fun_norm = @(T) tensor_norm(T);
        [nMPS_temp, Info] = LanczosExp(fun_prod, nMPS_temp, tau, 'innerprod', fun_inner, 'norm', fun_norm);
        nMPS2 = nMPS_temp;
        Kry_K2 = Info.K;
        Kry_K = max([Kry_K1, Kry_K2]);
    case 'L2R'
        MPS_temp = contract(MPS1, 2, MPS2, 1);
        %         MPS_temp = contract(MPS_temp, 1, H_envi_L, 2);
        %         MPS_temp = contract(MPS_temp, [1 5], H_local1, [4 1]);
        %         MPS_temp = contract(MPS_temp, [2 4], H_local2, [4 1]);
        %         MPS_temp = contract(MPS_temp, [1 4], H_envi_R, [2 3]);
        %         MPS_temp = permute(MPS_temp, [1 2 4 3]);
        
        fun_prod = @(T) twosite_prod(T, 'L2R', H_local1, H_local2, H_envi_R, H_envi_L);
        fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
        fun_norm = @(T) tensor_norm(T);
        [MPS_temp, Info] = LanczosExp(fun_prod, MPS_temp, -tau, 'innerprod', fun_inner, 'norm', fun_norm);
        Kry_K1 = Info.K;
        T_size = size(MPS_temp);
        MPS_mat = reshape(MPS_temp, T_size(1)*T_size(2), T_size(3)*T_size(4));
        [U, S, V, t_error, Se] = svdT(MPS_mat, 'Nkeep', D, 'epsilon', 1e-12);
        log_S_norm = log(norm(S));
        S = S/norm(S);
        S_size = size(S);
        nMPS1 = reshape(U, [], 2, S_size(1));
        nMPS1 = permute(nMPS1, [1 3 2]);
        nMPS2 = reshape(diag(S)*V, S_size(1), [], 2);
        
        H_envi = contract(H_envi_L, 2, nMPS1, 1);
        H_envi = contract(H_envi, [2 4], H_local1, [1 4]);
        H_envi = contract(H_envi, [1 4], conj(nMPS1), [1 3]);
        H_envi = permute(H_envi, [3 1 2]);
        nMPS_temp = nMPS2;
        %         nMPS_temp = contract(nMPS2, 1, H_envi, 2);
        %         nMPS_temp = contract(nMPS_temp, [2 4], H_local2, [4 1]);
        %         nMPS_temp = contract(nMPS_temp, [1 3], H_envi_R, [2 3]);
        %         nMPS_temp = permute(nMPS_temp, [1 3 2]);
        
        fun_prod = @(T) onesite_prod(T, H_local2, H_envi_R, H_envi);
        fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
        fun_norm = @(T) tensor_norm(T);
        [nMPS_temp, Info] = LanczosExp(fun_prod, nMPS_temp, tau, 'innerprod', fun_inner, 'norm', fun_norm);
        nMPS2 = nMPS_temp;
        Kry_K2 = Info.K;
        Kry_K = max([Kry_K1, Kry_K2]);
    case 'lL2R'
        MPS_temp = contract(MPS1, 2, MPS2, 1);
        %         MPS_temp = contract(MPS_temp, 3, H_local2, 3);
        %         MPS_temp = contract(MPS_temp, [2 3], H_local1, [4 2]);
        %         MPS_temp = contract(MPS_temp, [1 3], H_envi_L, [2 3]);
        %         MPS_temp = permute(MPS_temp, [3 2 1]);
        
        fun_prod = @(T) twosite_prod(T, 'lL2R', H_local1, H_local2, H_envi_R, H_envi_L);
        fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
        fun_norm = @(T) tensor_norm(T);
        [MPS_temp, Info] = LanczosExp(fun_prod, MPS_temp, -tau, 'innerprod', fun_inner, 'norm', fun_norm);
        Kry_K = Info.K;
        T_size = size(MPS_temp);
        MPS_mat = reshape(MPS_temp, T_size(1)*T_size(2), T_size(3));
        [U, S, V, t_error, Se] = svdT(MPS_mat, 'Nkeep', D, 'epsilon', 1e-12);
        log_S_norm = log(norm(S));
        S = S/norm(S);
        S_size = size(S);
        nMPS1 = reshape(U, [], 2, S_size(1));
        nMPS1 = permute(nMPS1, [1 3 2]);
        nMPS2 = reshape(diag(S)*V, S_size(1), 2);
        H_envi = [];
        
end
end