%============Testing coding================
clear; clc;
close all;

% addpath ../TensorFunction
% addpath ../SpinModel
% addpath ../stoMPS
addpath ./TDVP


Para.Geo.Lx = 3;
Para.Geo.Ly = 2;
Para.Geo.BCX = 'OBC';
Para.Geo.BCY = 'PBC';
Para.Model.J1xy = 1;
Para.Model.J1z = 1;
Para.Model.J2xy = 0;
Para.Model.J2z = 0;
Para.d = 2;
Para.L = Para.Geo.Lx * Para.Geo.Ly;
% [Intr] = IntrcMap_SLXXZ(Para);
[Intr] = IntrcMap_TLXXZ(Para);
[ H ] = AutomataInit( Para, Intr );


% beta = 0.5;
% tau = 0.125/16;
% [Para_Evo.beta_list] = Get_beta_expTau([2^-17, 0.5], [0.25/4, 10]);
[Para_Evo.beta_list] = Get_beta_expTau([2^-17, 1], [0.5/4, 10]);
Para_Evo.method = 'TDVP';
Para_Evo.D = 400;
Para_Evo.H_MPO = H;
Para_Evo.Lx = Para.Geo.Lx;
Para_Evo.Ly = Para.Geo.Ly;

Para_Evo.Phyval.Name = 'E_normal';
Para_Evo.Phyval.H = Para_Evo.H_MPO.A';
% Para_Evo.Phyval.Name = 'SpinCorMat';
% Para_Evo.Phyval.SpinType = 'Sz';

Para_Evo.SETTN.beta = Para_Evo.beta_list(1);
Para_Evo.SETTN.It_max = 20;
Para_Evo.SETTN.VariProd_step_max = 200;
Para_Evo.SETTN.MCrit = 100;
Para_Evo.SETTN.VariSum_step_max = 200;

Para_sto.Random_step = 2;  % Sampling number
Para_sto.D_0 = 1;          % Ini D_0
Para_sto.Rand_type = 'SON'; % Sampling space






% Filename = ['./Result/', 'Random_step=', num2str(Para_sto.Random_step),...
%     '_Sample_num_pertime=', num2str(Para_sto.D_0),...
%     '_D=', num2str(Para_Evo.D),'_tau=', num2str(tau),'_beta=', num2str(beta),...
%     '_N=',num2str(Para_Evo.Ly),'x', num2str(Para_Evo.Lx),'_Random_', Para_sto.Rand_type,'.mat'];





[Rslt_stoMPS, Info] = main_stoMPS(Para_sto, Para_Evo);


