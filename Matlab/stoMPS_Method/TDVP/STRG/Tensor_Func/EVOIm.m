function [ MPS_beta, anaHn, anaM ] = EVOIm( Para, MPS_0, Ham )
% input: Para, MPS_init, Ham
% output:MPS_beta = exp(-Para.beta/2 Ham)|MPS>, lognorm

beta = Para.beta;
MPS = MPS_0;
MPS_beta = MPS_0;

anaHn.err = zeros(Para.L-1, Para.It_max);
anaHn.entropy = zeros(Para.L-1, Para.It_max);

anaM.err = zeros(Para.L-1, Para.It_max);
anaM.entropy = zeros(Para.L-1, Para.It_max);
% keyboard;
Para.MCrit = 100;
for i = 1:1:Para.It_max
    %[MPS, lmp] = ProdMPS(MPS, Ham);     %H^i|j>
    [MPS.A, ns, anaHn1] = VariProdMPS(Para, MPS.A, Ham.A);     %H^i|j>
    
    anaHn.err(:,i) = anaHn1.err;
    anaHn.entropy(:,i) = anaHn1.entropy;
    
    MPS.lgnorm = MPS.lgnorm + log(ns) + Ham.lgnorm;
    
    sig = sign((-1)^i);
    fec1 = MPS_beta.lgnorm;
    fec2 = MPS.lgnorm + i * log(beta/2) - sum(log(1:1:i));
    % fprintf('%f, %f\n', lm/log(10), fec2/log(10));
    % fprintf('%d, fec1: %f, fec2: %f, delta: %f \n', i, fec1, fec2, abs(fec1-fec2)/log(10));
    if fec1 > fec2
        %[MPS_beta, lns] = MPSSum(MPS_beta, MPS, 1, sig * exp(fec2-fec1));
        [MPS_beta.A, ns, anaM1] = VariSumMPS(Para, MPS_beta.A, MPS.A, 1, sig * exp(fec2-fec1));
        MPS_beta.lgnorm = fec1 + log(ns);
    else
        %[MPS_beta, lns] = MPSSum(MPS_beta, MPS, exp(fec1-fec2), sig);
        [MPS_beta.A, ns, anaM1] = VariSumMPS(Para, MPS_beta.A, MPS.A, exp(fec1-fec2), sig);
        MPS_beta.lgnorm = fec2 + log(ns);
    end
    anaM.err(:,i) = anaM1.err;
    anaM.entropy(:,i) = anaM1.entropy;
%     [MPS_beta, lns] = MPSSum(MPS_beta, MPS, 1, 1);
%     lognorm = lognorm + lns;
    if (fec2-fec1)/log(10) < -16
        anaHn.err = anaHn.err(:, 1:1:i);
        anaHn.entropy = anaHn.entropy(:, 1:1:i);
        
        anaM.err = anaM.err(:, 1:1:i);
        anaM.entropy = anaM.entropy(:, 1:1:i);
%         load ED_res
%         MPSE = Cont(MPS_beta);
%         d1 = reshape(MPSE.data{1}, 2^Para.L, 1);
%         norm(d1 - ED_ans)
        break
    end
%     checkpsi( i, MPS, lm, MPS_beta, lognorm, Para )
    BX = bondcont(MPS_beta.A);
    fprintf('%d, %d, %.16f, %.16f\n', i, BX, fec2/log(10), (fec2-fec1)/log(10));
end

end

function BX = bondcont(MPS)
BX = 0;
for i = 1:1:length(MPS)
    if max(size(MPS{i})) > BX
        BX = max(size(MPS{i}));
    end
end
end