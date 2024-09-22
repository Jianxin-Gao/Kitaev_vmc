clear; clc;
close all;

Lx = 4;
Ly = 4;
site_num = Lx * Ly;
auto_cor_len = 20;

filename = '../../cmake-build-debug/one_point_functions';
file_id = fopen(filename,'rb');

row_data = fread(file_id,'double');


fclose(file_id);

sigma_row_data = row_data(1:site_num*2*3);
sigma_data = complex(sigma_row_data(1:2:end), sigma_row_data(2:2:end));

if length(row_data) == site_num * 9 + auto_cor_len * 2 % Number of processors > 1
    sigma_error_data = row_data(site_num*3*2+1:site_num*9);
    auto_cor_row_data = row_data(site_num*9+1:end);
    auto_cor_data = complex(auto_cor_row_data(1:2:end), auto_cor_row_data(2:2:end));
    
elseif length(row_data) == site_num * 3 * 2 + auto_cor_len * 2 % Number of processors = 1
    auto_cor_row_data = row_data(site_num*3*2+1:end);
    auto_cor_data = complex(auto_cor_row_data(1:2:end), auto_cor_row_data(2:2:end));
else
    error("Something wrong with number of data!");
end



