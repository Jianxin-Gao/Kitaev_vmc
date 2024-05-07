function [Rslt] = GetVal(Para_Phyval, MPS)
% Get observed values of physical quantities
% Para_Phyval.Name is a string that gives the name of the observed value
% e.g. 'SpinCorMat', 'Sz', 'G(beta/2)', 'Si_Sj'
Name = Para_Phyval.Name;
% try
    switch Name
        case 'SpinCorMat'
            % Get spin correlation matrix from MPS
            % Para_Phyval.SpinType, a string in {'Sx', 'Sy', 'Sz'} is required.
            SpinType = Para_Phyval.SpinType;
            Rslt = getCorReMat(MPS, SpinType);
            
        case 'Sz'
            % Get observed value of z-direction spin operator Sz
            % d = 2 only
            Rslt = getValOfSz(MPS);
            
        case 'E_normal'
            % Get energy from MPS (normalized)
            % Para_Phyval.H is required (H = H_MPO.A')
            H = Para_Phyval.H;
            Rslt = Energy(H, MPS);
    end
    
% catch
%     return;
% end