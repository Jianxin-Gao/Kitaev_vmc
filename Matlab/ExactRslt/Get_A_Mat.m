function A = Get_A_Mat(Model_Para, Intr)
% Return A_{j, k} = 2 * J_alpha{j, k} * u_{j, k}
% Model_Para.L = Lx =  Ly
% Model_Para.Jx = Jx, Model_Para.Jy = Jy;  Model_Para.Jz = Jz
% {j, k} is index, zigzag path

N = Model_Para.L^2;
A = zeros(N, N);
for i = 1:size(Intr, 1)
    j = Intr(i).Index1;
    k = Intr(i).Index2;
    J = Intr(i).CS;
    u = 1;
    A(j, k) = 2*J*u;
end
A = A - A';
end