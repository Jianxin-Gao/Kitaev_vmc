function [Rslt, Info] = main_stoMPS(Para_sto, Para_Evo)
maxNumCompThreads(4);

rng(sum(100*clock));
% rng(1);

%===================Stochastic parameter====================
Rand_type = Para_sto.Rand_type;
Random_step = Para_sto.Random_step;
D_0 = Para_sto.D_0;
%===================Evo parameter====================
Lx = Para_Evo.Lx;
Ly = Para_Evo.Ly;

switch Para_Evo.method
    case 'TEBD'
        Evo_func = @main_TEBD;
        K = length(Para_Evo.beta_list);
        flag_snake = 0;
    case 'TDVP'
        Evo_func = @main_TDVP_SETTN;
        K = length(Para_Evo.beta_list)-1; % -1 for SETTN
        flag_snake = 1;
end

%===================Rslt initialize====================
Rslt.Energy_Sum = zeros(1, K);
Rslt.Z_Sum = zeros(1, K);
Rslt.E_normal = zeros(Random_step, K);
Rslt.log_Z = zeros(Random_step, K);
Rslt.Z_list = zeros(Random_step, K);
PhyName = Para_Evo.Phyval.Name; SumName = [PhyName, '_Sum'];
Rslt.(PhyName) = cell(Random_step, K);
Rslt.(SumName) = num2cell(zeros(1, K));
Info.Se = zeros(Random_step, K);
Info.Kry_K = zeros(Random_step, K);
Info.trun_err = zeros(Random_step, K);
%==================begin sampling=======================
for j = 1:1:Random_step
    fprintf([repmat('=', 1, 34),'Random step: %d/%d', repmat('=', 1, 34), '\n'], j, Random_step)
    
    if flag_snake
        MPS = Initialize(Lx * Ly, Rand_type, D_0, 2);
    else
        MPS = Initialize(Lx, Rand_type, D_0, 2^Ly);
    end
    
    A_cell = Canonicalization(MPS);
    
    [Rslt_step, Info_step] =  Evo_func(Para_Evo, A_cell);
    %     Rslt.MPS(j, :) = Rslt_step.MPS;
    beta_list = Rslt_step.beta_list;
    Rslt.log_Z(j, :) = Rslt_step.log_Z;
    Rslt.Z_list(j, :) = exp(Rslt_step.log_Z);
    Rslt.Z_Sum = Rslt.Z_Sum + Rslt.Z_list(j, :);
    
    Rslt.(PhyName)(j, :) = Rslt_step.(PhyName);
    temp = arrayfun(@(x)Rslt_step.(PhyName){x}.*Rslt.Z_list(j, x), 1:K, 'UniformOutput', false);
    Rslt.(SumName) = arrayfun(@(y)(Rslt.(SumName){y}+temp{y}), 1:K, 'UniformOutput', false);
    
    Info.Se(j, :) = Info_step.Se;
    Info.Kry_K(j, :) = Info_step.Kry_K;
    Info.trun_err(j, :) = Info_step.trun_err;
end

[RsltName] = getName(PhyName);
%[eta] = Rand_eff(Z_list, Z_Sum);
Rslt.Z_Rslt = Rslt.Z_Sum./Random_step;
Rslt.(RsltName) = arrayfun(@(z)(Rslt.(SumName){z}./Rslt.Z_Sum(z)), 1:K, 'UniformOutput', false);
Rslt.beta_list = beta_list;
end

function [RsltName] = getName(PhyName)
switch PhyName
    case 'E_normal'
        RsltName = 'Energy_Rslt';
    otherwise
        RsltName = [PhyName, '_Rslt'];
end
end