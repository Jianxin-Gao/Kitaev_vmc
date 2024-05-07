%===============Test coding====================
clear; clc;
% profile on
addpath ../TensorFunction
addpath ../SpinModel
addpath ../stoMPS

Para.D = 50;
n = 5; D_0 = 1;
% Para.H_MPO = Automata(n, 1);

Para.H_MPO.A = Automata(n, 1);
Para.H_MPO.lgnorm = 0;
Para.J = 1;
Para.spin_model = 'Heisenberg';
iniMPS = getMPS_sphere(n, D_0, 2);

beta_begin = 1/2^31;
beta_list = [beta_begin];
for i = 1:round(log2(1/beta_begin))
    beta_list = [beta_list 2*beta_list(end)];
end
for j = 1:21
    beta_list = [beta_list 2*j+1];
end

Para.beta_list = beta_list;

Para.SETTN.beta = Para.beta_list(1);
Para.SETTN.It_max = 20;
Para.SETTN.VariProd_step_max = 200;
Para.SETTN.MCrit = 10000;
Para.SETTN.VariSum_step_max = 200;
% 
Para.Phyval.Name = 'E_normal';
Para.Phyval.H = Para.H_MPO.A';

[Rslt_tdvp] = main_TDVP_SETTN(Para, iniMPS);
% [Rslt_tdvp] = main_TDVP(Para, iniMPS);
% Rslt_tdvp.E = cell2mat(Rslt_tdvp.E_normal);

[Rslt_ED] = ED_MPS(Para, iniMPS);

beta_list = beta_list(2:end);
Rslt_ED.E = Rslt_ED.E(2:end);
plot(beta_list, abs((cell2mat(Rslt_tdvp.E) - Rslt_ED.E)./Rslt_ED.E), '-o')
% profile viewer