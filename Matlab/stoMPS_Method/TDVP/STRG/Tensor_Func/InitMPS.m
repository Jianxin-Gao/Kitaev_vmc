function [ MPS ] = InitMPS( Para, str, dir )
% Initialize the MPS
%   |     |     |     |
%   ^     ^     ^     ^
%   |     |     |     |
%   o--<--o--<--o--<--o

MPS.A = cell(Para.L, 1);
MPS.lgnorm = 0;
for i = 1:1:Para.L
    switch i
        case {1, Para.L}
            MPS.A{i} = LocS(str(i), dir);
            MPS.A{i} = reshape(MPS.A{i}, [1,2]);
        otherwise
            MPS.A{i} = LocS(str(i), dir);
            MPS.A{i} = reshape(MPS.A{i}, [1,1,2]);
    end
end

end

function vec = LocS(str, dir)
switch str
    case '1'
        vec = [1;0];
    case '-1'
        vec = [0;1];
end
if strcmp(dir, 'x')
    vec = 1/sqrt(2) * [1,1;1, -1] * vec;
end
end