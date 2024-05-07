clear; clc;

addpath /Users/jxgao/Library/CloudStorage/OneDrive-个人/jxGao_project/2023-stoMPS/Code/stoMPS_2D
addpath /Users/jxgao/Library/CloudStorage/OneDrive-个人/jxGao_project/2023-stoMPS/Code/SpinModel
addpath /Users/jxgao/Library/CloudStorage/OneDrive-个人/jxGao_project/2023-stoMPS/Code/TensorFunction
addpath /Users/jxgao/Library/CloudStorage/OneDrive-个人/jxGao_project/2023-stoMPS/Code/TDVP

beta_max = 0.5;
[~, beta_list ] = GetStepFunc(2^-15, [beta_max, 0.125], []);

Para.beta = 2*beta_list(1);
load('IniMPS-chi=20-d=2-n=32.mat');
MPS_0.A = iniMPS;
MPS_0.lgnorm = 0;
Para.L = 32;
Para.It_max = 1;
Para.Geo.Lx = 8;
Para.Geo.Ly = 4;
Para.Geo.BCY = 'PBC';
Para.Model.J1xy = 1;
Para.Model.J1z = 1;
Para.d = 2;
[ Intr ] = IntrcMap_SLXXZ( Para );
[ H ] = AutomataInit( Para, Intr );
Ham = H;
Para.VariProd_step_max = 200;
Para.MCrit = 300;
Para.VariSum_step_max = 200;


Para_sto.Random_step = 1;
Para_Evo.beta_list = beta_list;
Para_Evo.method = 'TDVP';
Para_Evo.D = 300;
Para_Evo.H_MPO = Ham.A';
[Rslt] = main_stoMPS_fixinimps(Para_sto, Para_Evo, iniMPS);



% [ MPS_beta, anaHn, anaM ] = EVOIm( Para, MPS_0, Ham );
% 
% Para.D = 300;
% Para.H_MPO = Ham.A;
% MPS_beta.A{1} = MPS_beta.A{1} * exp(MPS_beta.lgnorm);
% iniMPS = MPS_beta.A;
% 
% Para_sto.Random_step = 1;
% Para_Evo.beta_list = beta_list(2:end);
% Para_Evo.method = 'TDVP';
% Para_Evo.D = 300;
% Para_Evo.H_MPO = Ham.A';
% [Rslt] = main_stoMPS_fixinimps(Para_sto, Para_Evo, iniMPS);
