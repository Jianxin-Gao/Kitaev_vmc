function [Rslt] = GetED_Ham(beta_lis, Ham_mat)
K = length(beta_lis);
%====================get Energy & C=======================
Rslt.characteristic_value = eig(Ham_mat);
Rslt.beta = beta_lis;
Rslt.Z = zeros(1, K);
Rslt.E_Z = zeros(1, K);

Rslt.C = zeros(1, K);
for i = 1:K
    beta = beta_lis(i);
    Rslt.Z(i) = sum(exp(-beta*Rslt.characteristic_value));
    Rslt.E_Z(i) = sum(exp(-beta*Rslt.characteristic_value).*Rslt.characteristic_value);
    Rslt.C(i) = -beta^2 * (-sum(exp(-beta*Rslt.characteristic_value).*Rslt.characteristic_value.^2)/Rslt.Z(i) + (Rslt.E_Z(i)/Rslt.Z(i))^2);
end
Rslt.E = Rslt.E_Z./Rslt.Z;
end


