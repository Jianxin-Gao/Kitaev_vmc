function [ ] = Main( Para, time, init_state )
%MAIN Summary of this function goes here
%   Detailed explanation goes here
rng('shuffle', 'twister')
str = zeros(Para.L, 1);

for i = 1:1:Para.L
    switch init_state
        case 'FM'
            str(i) = 1;
        case 'AF'
            str(i) = (-1)^i;
    end
end
str_init = str;
MPS = InitMPS(Para, string(str), 'o');
Enlist = zeros(Para.Step_max, 1);
aveEnlist = zeros(Para.Step_max, 1);
errorbarlist = zeros(Para.Step_max, 1);
tic
[Ham, ~, ~] = InitHam(Para); 
toc
fprintf('=============================================================== \n');
for i = 1:1:Para.Step_max
    if mod(i, 2) == 1
        dir = 'x';
    else
        dir = 'z';
    end
    tic
    [ MPS_e, anaHn, anaM ] = EVOIm( Para, MPS, Ham );
    toc
    if Para.Save_EE == 1
        save(['anaRslt/EE_Rslt/', Para.EE_name, num2str(i), '.mat'], 'anaHn', 'anaM', 'Para');
    end
%     keyboard;
%     if i == 1
%         overlap = zeros(10,2);
%         METS_init = MPS_e;
%         overlap(i, 1) = i;
%         overlap(i, 2) = overlap_check(METS_init, MPS_e);
%         fprintf('%5d | %.8f \n', i, overlap(i, 2));
%     else
%         overlap(i, 1) = i;
%         overlap(i, 2) = overlap_check(METS_init, MPS_e);
%         fprintf('%5d | %.8f \n', i, overlap(i, 2));
%     end
%     if i == 2
%         keyboard;
%     end
    Enlist(i) = real(Observe(MPS_e, Ham, Para));
%     keyboard;
    str = cell(50, 1);
    for it = 1:1:50
    str{it} = Sampling(MPS_e.A, dir);
    
    end
    save('str.mat', 'str');
    keyboard;
    MPS = InitMPS(Para, string(str), 'dir');
    Bx = bondcont(MPS_e.A);
    aveEnlist(i) = sum(Enlist(1:1:i))/i;
    aveval = sum(aveEnlist(1:1:i))/i;
    errorbarlist(i) = sum((aveEnlist(1:1:i)-aveval).^2)/i;
    fprintf('|Number of samples: %6d| Bond: %3d| <En>: %.8f| Error bar: %8f| \n', i, Bx, aveEnlist(i), errorbarlist(i));
    if mod(i, 500) == 0
        save(['SamRslt/Enlist-beta=', num2str(Para.beta),'-Ns=', num2str(i), '-', time, '.mat'], 'Enlist', 'aveEnlist', 'errorbarlist', 'str_init', 'str');
    end
    if i == Para.Step_max
        save(['anaRslt/EE_Rslt/', Para.EE_name, 'end', '.mat'], 'str');
    end
end
fprintf('=============================================================== \n');
end


function BX = bondcont(MPS)
BX = 0;
for i = 1:1:length(MPS)
    if max(size(MPS{i})) > BX
        BX = max(size(MPS{i}));
    end
end
end
