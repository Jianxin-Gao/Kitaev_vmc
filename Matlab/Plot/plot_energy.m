clear; clc;
close all;

Lx = 4;
Ly = 4;


filename = '../../cmake-build-debug/energy_statistics';
file_id = fopen(filename,'rb');

row_data = fread(file_id,'double');

fclose(file_id);

% Get energy
energy = complex(row_data(1), row_data(2));

% Get energy error
energy_err = row_data(3);

% Get bond enenrgy
Num_horizontal_bond_energy = (Lx - 1) * Ly;
Num_vertical_bond_energy = (Ly - 1) * Lx;
horizontal_bond_energy_row_data = row_data(4:4+Num_horizontal_bond_energy*2-1);
horizontal_bond_energy = complex(horizontal_bond_energy_row_data(1:2:end), horizontal_bond_energy_row_data(2:2:end));
vertical_bond_energy_row_data = row_data(4+Num_horizontal_bond_energy*2:4+Num_horizontal_bond_energy*2+Num_vertical_bond_energy*2-1);
vertical_bond_energy = complex(vertical_bond_energy_row_data(1:2:end), vertical_bond_energy_row_data(2:2:end));

