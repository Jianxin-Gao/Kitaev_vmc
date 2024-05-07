load('str.mat')
[ Para ] = GetPara( 5 );
[ MPS ] = InitMPS( Para, str{1}, 'o' );
for i = 2:1:49
A = InitMPS( Para, str{i}, 'o' );
%keyboard;
MPS.A = DrctSumMPS(MPS.A, A.A);
end
A = InitMPS( Para, str{end}, 'o' );
[ C, Ns2, ana ] = VariSumMPS(Para, MPS.A, A.A);