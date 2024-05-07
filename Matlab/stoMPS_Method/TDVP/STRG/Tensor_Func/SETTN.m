function [ rho, Op, anaHn, anaM ] = SETTN( Para, H, Id, Op )
% function [ rho, Op ] = SETTN( Para, H, Id, Op )
% Series-Expansion Thermal Tensor Network algorithm to get 
% the density matrix rho at inverse temperature Para.tau.
% rho.A is the tensor that make up the MPO of density martrix.
% rho.lgnorm is the log norm of density martrix.
% Yuan Gao@buaa 2020.12.07
% mail: 17231064@buaa.edu.cn

tau = Para.beta/2;

rho = Id;
Hn = H;

anaHn.err = zeros(Para.L-1, Para.It_max);
anaHn.entropy = zeros(Para.L-1, Para.It_max);

anaM.err = zeros(Para.L-1, Para.It_max);
anaM.entropy = zeros(Para.L-1, Para.It_max);

for i = 1:1:Para.It_max
    % //rho = rho + (-tau)^i/i! H^i
    sig = sign((-1)^i);
    fec1 = rho.lgnorm;
    fec2 = Hn.lgnorm + i * log(tau) - sum(log(1:1:i));
    
    if fec1 > fec2
        [rho.A, nm, anaM1] = VariSumMPO(Para, rho.A, Hn.A, [1, sig * exp(fec2-fec1)]);
        rho.lgnorm = fec1 + log(nm);
    else
        [rho.A, nm, anaM1] = VariSumMPO(Para, rho.A, Hn.A, [exp(fec1-fec2), sig]);
        rho.lgnorm = fec2 + log(nm);
    end
    fprintf('%d, %.5f, %.5f\n', i, fec1/log(10), fec2/log(10));
    if (fec2-fec1)/log(10) < -16
        anaHn.err = anaHn.err(:, 1:1:i);
        anaHn.entropy = anaHn.entropy(:, 1:1:i);
        
        anaM.err = anaM.err(:, 1:1:i);
        anaM.entropy = anaM.entropy(:, 1:1:i);
        break;
    end
    anaM.err(:,i) = anaM1.err;
    anaM.entropy(:,i) = anaM1.entropy;
    % //Hn = H * Hn (H^(i+1))
    [Hn.A, nm, anaHn1] = VariProdMPO(Para, Hn.A, H.A);
    anaHn.err(:,i) = anaHn1.err;
    anaHn.entropy(:,i) = anaHn1.entropy;
    Hn.lgnorm = H.lgnorm + Hn.lgnorm + log(nm);
end
end

