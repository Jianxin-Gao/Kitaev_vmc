% addpath tensor_func
% addpath ED
% Model.J = 1;
% Model.alpha = 1;
% Para = GetPara(1);
% Para.Model = Model;
% [ ED_Rslt ] = GetEDRslt( Para, 1 );
% save('ED_res.mat', 'ED_Rslt')
% Run()
% addpath Tensor_Func
% beta = 5;
% load('MPS.mat')
% load('Ham.mat')
% v = MPS2V(MPS.A);
% res_RG = MPS2V(MPS_e.A) * expm(MPS_e.lgnorm);
% norm(res_RG - expm(-beta/2 * H) * v)/norm(expm(-beta/2 * H) * v)
% transp(expm(-beta/2 * H) * v) * H * expm(-beta/2 * H) * v/(transp(expm(-beta/2 * H) * v) * expm(-beta/2 * H) * v)
% trace(expm(-beta * H) * H)/trace(expm(-beta * H))