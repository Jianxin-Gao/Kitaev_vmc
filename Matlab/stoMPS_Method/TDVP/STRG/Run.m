function [ ] = Run()

addpath Tensor_Func

maxNumCompThreads(5);

beta = 5; 
TStr = datestr(now,'YYYYmmDD_HHMMSS');
Para = GetPara(beta);
Main(Para, TStr, 'AF');

end