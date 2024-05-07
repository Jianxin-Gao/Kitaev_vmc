function [A_cell, theta_cell] = Initialize(n, type, varargin)
%varargin{1} for D_0, {} for D_0 = 1;
%varargin{2} for d, {} for d = 2;
switch type
    case 'FM_groundState'
        A = ones(1,1,2);
        A(1,1,2) = 0;
        A_bond = ones(1,2);
        A_bond(1,2) = 0;
        A_cell = {A_bond};
        for i = 2:(n-1)
            A_cell{i} = A;
        end
        A_cell{n} = A_bond;
        theta_cell = {};
    case 'IntrinsicSxAll'
        [Sx, ~, ~, ~] = SpinOp(2);
        [Vec, ~] = eig(Sx);
        
        A = ones(1,1,2);
        A(1, 1, 1) = Vec(1, 1);
        A(1,1,2) = Vec(2, 1);
        A_bond = ones(1,2);
        A_bond(1, 1) = Vec(1, 1);
        A_bond(1,2) = Vec(2, 1);
        A_cell = {A_bond};
        for i = 2:(n-1)
            A_cell{i} = A;
        end
        A_cell{n} = A_bond;
        theta_cell = {};
    case  'Random_IntrinsicSx'
        [Sx, ~, ~, ~] = SpinOp(2);
        [Vec, ~] = eig(Sx);
        A_bond_L = ones(1, 2);
        random = randi(2)-1;
        if random == 1
            A_bond_L(1, 1) = Vec(1, 1);
            A_bond_L(1, 2) = Vec(2, 1);
        else
            A_bond_L(1, 1) = Vec(1, 2);
            A_bond_L(1, 2) = Vec(2, 2);
            
        end
        A_cell = {A_bond_L};
        for i = 2:(n-1)
            A = ones(1, 1, 2);
            random = randi(2)-1;
            if random == 1
                A(1, 1, 1) = Vec(1, 1);
                A(1, 1, 2) = Vec(2, 1);
            else
                A(1, 1, 1) = Vec(1, 2);
                A(1, 1, 2) = Vec(2, 2);
            end
            A_cell{i} = A;
        end
        A_bond_R = ones(1, 2);
        random = randi(2)-1;
        if random == 1
            A_bond_R(1, 1) = Vec(1, 1);
            A_bond_R(1, 2) = Vec(2, 1);
        else
            A_bond_R(1, 1) = Vec(1, 2);
            A_bond_R(1, 2) = Vec(2, 2);
        end
        A_cell{n} = A_bond_R;
        theta_cell = {};
        
    case  'SO2'
        if isempty(varargin)
            theta = rand*2*pi;
            theta_cell = {theta};
            A_bond_L = ones(1, 2);
            A_bond_L(1, 1) = cos(theta);
            A_bond_L(1, 2) = sin(theta);
            A_cell = {A_bond_L};
            for i = 2:(n-1)
                theta = rand*2*pi;
                theta_cell{i} = theta;
                A = ones(1, 1, 2);
                A(1, 1, 1) = cos(theta);
                A(1, 1, 2) = sin(theta);
                A_cell{i} = A;
            end
            theta = rand*2*pi;
            theta_cell{n} = theta;
            A_bond_R = ones(1, 2);
            A_bond_R(1, 1) = cos(theta);
            A_bond_R(1, 2) = sin(theta);
            A_cell{n} = A_bond_R;
        elseif varargin{1} == 1
            [A_cell, theta_cell] = Initialize(n, type);
        else
            theta_cell = cell(varargin{1}, n);
            [A_cell, theta_cell(1, :)] = Initialize(n, type);
            for l = 1:1:varargin{1}-1
                [tempA, theta_cell(l+1, :)] = Initialize(n, type);
                A_cell = direct_sum(A_cell, tempA);
            end
        end
        
    case  'Z2'
        if isempty(varargin)
            A_bond_L = ones(1, 2);
            random = randi(2)-1;
            if random == 1
                A_bond_L(1, 1) = 0;
            else
                A_bond_L(1, 2) = 0;
            end
            A_cell = {A_bond_L};
            for i = 2:(n-1)
                A = ones(1, 1, 2);
                random = randi(2)-1;
                if random == 1
                    A(1, 1, 1) = 0;
                else
                    A(1, 1, 2) = 0;
                end
                A_cell{i} = A;
            end
            A_bond_R = ones(1, 2);
            random = randi(2)-1;
            if random == 1
                A_bond_R(1, 1) = 0;
            else
                A_bond_R(1, 2) = 0;
            end
            A_cell{n} = A_bond_R;
            
        elseif varargin{1} == 1
            A_cell = Initialize(n, type);
        else
            A_cell = Initialize(n, type);
            for l = 1:1:varargin{1}-1
                A_cell = direct_sum(A_cell, Initialize(n, type));
            end
        end
        theta_cell = {};
        
    case 'SON'
        if isempty(varargin)
            [A_cell, theta_cell] = Initialize(n, 'SO2');
        elseif ~isempty(varargin) && varargin{1} ~= 1 && length(varargin) == 2
            A_cell = getMPS_sphere(n, varargin{1}, varargin{2});
            theta_cell = {};
        elseif ~isempty(varargin) && varargin{1} ~= 1 && length(varargin) == 1
            A_cell = getMPS_sphere(n, varargin{1}, 2);
            theta_cell = {};
        elseif varargin{1} == 1
            [A_cell, theta_cell] = Initialize(n, 'SO2');
        end
        
end

end
