function [Rslt] = kron_continuously(matrix_cell)
Rslt = kron(matrix_cell{1}, matrix_cell{2});
matrix_cell_size = size(matrix_cell);
n = matrix_cell_size(2);
for i = 2:1:(n-1)
    Rslt = kron(Rslt, matrix_cell{i+1});
end
end