function Ham_mat = GetHam_IntrMap(IntrMap, n)
Ham_mat = Get_Ham_site(IntrMap(1), n);
for j = 2:size(IntrMap, 1)
    tic
    Ham_mat = Ham_mat + Get_Ham_site(IntrMap(j), n);
    fprintf('%d/%d', j, size(IntrMap, 1))
    toc
end
end

function [mat_cell] = Initial_Id(n, Id)
%输出n个全Id组成的cell
mat_cell = num2cell(ones(1, n));
for i = 1:n
    mat_cell{i} = Id;
end
end

function [H_site] = Get_Ham_site(Intr_site, n)
d = 2;
[Sx, Sy, Sz, Id] = SpinOp(d);

temp = Initial_Id(n, Id);
switch Intr_site.JmpOut_type
    case 'Sx'
        temp{Intr_site.JmpOut} = Sx;
    case 'Sy'
        temp{Intr_site.JmpOut} = Sy;
    case 'Sz'
        temp{Intr_site.JmpOut} = Sz;
end

switch Intr_site.JmpIn_type
    case 'Sx'
        temp{Intr_site.JmpIn} = Sx;
    case 'Sy'
        temp{Intr_site.JmpIn} = Sy;
    case 'Sz'
        temp{Intr_site.JmpIn} = Sz;
end

H_site = Intr_site.CS * kron_continuously(temp);

end