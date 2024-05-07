clear; clc;
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -1j; 1j, 0];
sigma_z = [1, 0; 0, -1];

spin_up = [1; 0];
spin_down = [0; 1];

S_plus = 1/2 * (sigma_x + 1j * sigma_y);
S_minus = 1/2 * (sigma_x - 1j * sigma_y);

S_plus * spin_down

sigma_x * spin_down

yy = kron(sigma_y, sigma_y);
yy * kron(spin_down, spin_down)

sigma_z * spin_up
zz = kron(sigma_z, sigma_z)
zz * kron(spin_up, spin_up)

xx = kron(sigma_x,sigma_x)
xx * kron(spin_up, spin_down)
xx * kron(spin_up, spin_up)