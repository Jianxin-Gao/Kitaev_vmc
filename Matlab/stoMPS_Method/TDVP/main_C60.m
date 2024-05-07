clear; clc;
Para.n = 90;
Para.D = 32;
beta_begin = 1/2^31;
beta_list = [beta_begin];
for i = 1:round(log2(1/beta_begin))
    beta_list = [beta_list 2*beta_list(end)];
end
for j = 1:21
    beta_list = [beta_list 2*j+1];
end
Para.beta_list = beta_list;

load('./Rslt_C60/H_theta=0.7500pi.mat');
Para.H_MPO = H;

[Rslt] = main_TDVP(Para);
[MPS_METTS] = Sample(Rslt.MPS);
[Sz_val] = getValOfSz(MPS_METTS);
[E] = Energy(H_MPO, MPS_METTS);

%save('./Rslt_C60/C60_groundState_theta=0.7500pi.mat', 'Sz_val', 'E')