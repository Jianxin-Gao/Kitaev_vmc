% addpath tensor_func

Para.VariSum_step_max = 1000;
Para.VariProd_step_max = 1000;
Para.MCrit = 200;
ep = 1e-8;



A = cell(10,1);
B = cell(10,1);

A{1} = rand(4,2) + 1i * rand(4,2);
A{end} = rand(4,2) + 1i * rand(4,2);
B{1} = rand(4,2,2) + 1i * rand(4,2,2);
B{end} = rand(4,2,2) + 1i * rand(4,2,2);
for i = 2:1:(length(A)-1)
    A{i} = rand(4,4,2) + 1i * rand(4,4,2);
    B{i} = rand(4,4,2,2) + 1i * rand(4,4,2,2);
end

% tic
[C, Ns] = VariProdMPS(Para, A,B);
% toc
Cp = ProdMPS(A,B);

T1 = MPS2V(C);
T2 = MPS2V(Cp);
norm(T1*Ns-T2)/4^length(A)/Ns