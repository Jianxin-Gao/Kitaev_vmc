addpath Tensor_Func
addpath Model
addpath ED
Model.Delta = 1;
Model.alpha = 1;
beta = 5;
Para = GetPara(beta);
Para.Model = Model;
tic
[H, Id, Op] = InitHam(Para);
[rho_init, Op, anaHn, anaRho] = SETTN(Para, H, Id, Op);
save('anaRslt/SETTNRES.mat', 'anaHn', 'anaRho', 'Para')
load('Ham.mat')
rho = MPO2H(rho_init.A) * exp(rho_init.lgnorm);
norm(rho-expm(-beta/2*H))/norm(expm(-beta/2*H))