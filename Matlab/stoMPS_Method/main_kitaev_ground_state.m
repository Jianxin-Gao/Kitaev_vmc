clear; clc;
addpath ./TDVP/

Para.Model.J1x = 1;
Para.Model.J1y = 1;
Para.Model.J1z = 1;
Para.L = 4*4;
Para.d = 2;

[ Intr ] = IntrcMap_Kitaev(Para);
[ H ] = AutomataInit( Para, Intr );

[Para_Evo.beta_list] = Get_beta_expTau([2^-17, 1], [0.5/2, 100]);
Para_Evo.method = 'TDVP';
Para_Evo.D = 50;
Para_Evo.H_MPO = H;
Para_Evo.Lx = 4;
Para_Evo.Ly = 4;

Para_Evo.Phyval.Name = 'E_normal';
Para_Evo.Phyval.H = Para_Evo.H_MPO.A';
% Para_Evo.Phyval.Name = 'SpinCorMat';
% Para_Evo.Phyval.SpinType = 'Sz';

Para_Evo.SETTN.beta = Para_Evo.beta_list(1);
Para_Evo.SETTN.It_max = 20;
Para_Evo.SETTN.VariProd_step_max = 200;
Para_Evo.SETTN.MCrit = 100;
Para_Evo.SETTN.VariSum_step_max = 200;

Para_sto.Random_step = 1;  % Sampling number
Para_sto.D_0 = 1;          % Ini D_0
Para_sto.Rand_type = 'SON'; % Sampling space

[Rslt_stoMPS, Info] = main_stoMPS(Para_sto, Para_Evo);


