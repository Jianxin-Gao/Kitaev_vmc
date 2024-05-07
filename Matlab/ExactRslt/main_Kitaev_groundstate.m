clear;
Model_Para.L = 4;
Model_Para.Jx = 1;
Model_Para.Jy = 1;
Model_Para.Jz = 0.01;

[ Intr ] = IntrcMap_Kitaev(Model_Para);
A = Get_A_Mat(Model_Para, Intr);

E = -sum(abs(eig(A)))/2/2;
En = E/Model_Para.L^2;