function [ ED_Rslt ] = GetEDRslt( Para, beta_list )
InMap = IntrcMap(Para);
[H, M] = ED_Hamiltonian(Para, InMap);
save('Ham.mat', 'H');
[~, Cm, LnZ, En] = ED_Cm(H, beta_list);
ED_Rslt.Cm = Cm/Para.L;
ED_Rslt.LnZ = LnZ/Para.L;
ED_Rslt.En = En/Para.L;
if norm(Para.Field.h) ~= 0
    [~, ~, M] = ED_chi(H, M, norm(Para.Field.h), beta_list);
    ED_Rslt.M = real(M)/Para.L;
end
end

