function [E] = Energy(H_MPO, A_cell)
A_cell_size = size(A_cell);
n = A_cell_size(2);
H_cell = H_MPO;

H_cont = contract(H_cell{1}, 3, A_cell{1}, 2);
H_cont = contract(H_cont, 2, conj(A_cell{1}), 2);
H_cont = permute(H_cont, [3 2 1]);
for i = 2:1:(n-1)
    H_cont = contract(H_cont, 2, A_cell{i}, 1);
    H_cont = contract(H_cont, [2 4], H_cell{i}, [1 4]);
    H_cont = contract(H_cont, [1 4], conj(A_cell{i}), [1 3]);
    H_cont = permute(H_cont, [3 1 2]);
end
H_cont = contract(H_cont, 2, A_cell{n}, 1);
H_cont = contract(H_cont, [2 3], H_cell{n}, [1 3]);
H_cont = contract(H_cont, [1 2], conj(A_cell{n}), [1 2]);
E = H_cont;
end

