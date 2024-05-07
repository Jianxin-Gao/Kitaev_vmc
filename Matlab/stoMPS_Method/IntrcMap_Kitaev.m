function [ Intr ] = IntrcMap_Kitaev(Para)
% for 4x4 Kitaev

int_num = 18;
int_cell = cell(int_num, 5);
count = 1;
% ========== 1 ============
int_cell{count, 1} = 1;
int_cell{count, 2} = 5;
int_cell{count, 3} = 'Sx';
int_cell{count, 4} = 'Sx';
int_cell{count, 5} = Para.Model.J1x;
count = count + 1;
% ========== 2 ============
int_cell{count, 1} = 9;
int_cell{count, 2} = 13;
int_cell{count, 3} = 'Sx';
int_cell{count, 4} = 'Sx';
int_cell{count, 5} = Para.Model.J1x;
count = count + 1;
% ========== 3 ============
int_cell{count, 1} = 6;
int_cell{count, 2} = 10;
int_cell{count, 3} = 'Sx';
int_cell{count, 4} = 'Sx';
int_cell{count, 5} = Para.Model.J1x;
count = count + 1;
% ========== 4 ============
int_cell{count, 1} = 3;
int_cell{count, 2} = 7;
int_cell{count, 3} = 'Sx';
int_cell{count, 4} = 'Sx';
int_cell{count, 5} = Para.Model.J1x;
count = count + 1;
% ========== 5 ============
int_cell{count, 1} = 11;
int_cell{count, 2} = 15;
int_cell{count, 3} = 'Sx';
int_cell{count, 4} = 'Sx';
int_cell{count, 5} = Para.Model.J1x;
count = count + 1;
% ========== 6 ============
int_cell{count, 1} = 8;
int_cell{count, 2} = 12;
int_cell{count, 3} = 'Sx';
int_cell{count, 4} = 'Sx';
int_cell{count, 5} = Para.Model.J1x;
count = count + 1;

% ========== 7 ============
int_cell{count, 1} = 2;
int_cell{count, 2} = 3;
int_cell{count, 3} = 'Sy';
int_cell{count, 4} = 'Sy';
int_cell{count, 5} = Para.Model.J1y;
count = count + 1;
% ========== 8 ============
int_cell{count, 1} = 5;
int_cell{count, 2} = 6;
int_cell{count, 3} = 'Sy';
int_cell{count, 4} = 'Sy';
int_cell{count, 5} = Para.Model.J1y;
count = count + 1;
% ========== 9 ============
int_cell{count, 1} = 7;
int_cell{count, 2} = 8;
int_cell{count, 3} = 'Sy';
int_cell{count, 4} = 'Sy';
int_cell{count, 5} = Para.Model.J1y;
count = count + 1;
% ========== 10 ============
int_cell{count, 1} = 10;
int_cell{count, 2} = 11;
int_cell{count, 3} = 'Sy';
int_cell{count, 4} = 'Sy';
int_cell{count, 5} = Para.Model.J1y;
count = count + 1;
% ========== 11 ============
int_cell{count, 1} = 13;
int_cell{count, 2} = 14;
int_cell{count, 3} = 'Sy';
int_cell{count, 4} = 'Sy';
int_cell{count, 5} = Para.Model.J1y;
count = count + 1;
% ========== 12 ============
int_cell{count, 1} = 15;
int_cell{count, 2} = 16;
int_cell{count, 3} = 'Sy';
int_cell{count, 4} = 'Sy';
int_cell{count, 5} = Para.Model.J1y;
count = count + 1;

% ========== 13 ============
int_cell{count, 1} = 1;
int_cell{count, 2} = 2;
int_cell{count, 3} = 'Sz';
int_cell{count, 4} = 'Sz';
int_cell{count, 5} = Para.Model.J1z;
count = count + 1;
% ========== 14 ============
int_cell{count, 1} = 3;
int_cell{count, 2} = 4;
int_cell{count, 3} = 'Sz';
int_cell{count, 4} = 'Sz';
int_cell{count, 5} = Para.Model.J1z;
count = count + 1;
% ========== 15 ============
int_cell{count, 1} = 6;
int_cell{count, 2} = 7;
int_cell{count, 3} = 'Sz';
int_cell{count, 4} = 'Sz';
int_cell{count, 5} = Para.Model.J1z;
count = count + 1;
% ========== 16 ============
int_cell{count, 1} = 9;
int_cell{count, 2} = 10;
int_cell{count, 3} = 'Sz';
int_cell{count, 4} = 'Sz';
int_cell{count, 5} = Para.Model.J1z;
count = count + 1;
% ========== 17 ============
int_cell{count, 1} = 11;
int_cell{count, 2} = 12;
int_cell{count, 3} = 'Sz';
int_cell{count, 4} = 'Sz';
int_cell{count, 5} = Para.Model.J1z;
count = count + 1;
% ========== 18 ============
int_cell{count, 1} = 14;
int_cell{count, 2} = 15;
int_cell{count, 3} = 'Sz';
int_cell{count, 4} = 'Sz';
int_cell{count, 5} = Para.Model.J1z;
count = count + 1;

Intr = struct('JmpOut', int_cell(:,1), 'JmpIn', int_cell(:,2), ...
    'JmpOut_type', int_cell(:,3), 'JmpIn_type', int_cell(:,4), 'CS', int_cell(:,5));

for i = 1:1:length(Intr)
    if Intr(i).JmpOut > Intr(i).JmpIn
        T = Intr(i).JmpIn;
        Intr(i).JmpIn = Intr(i).JmpOut;
        Intr(i).JmpOut = T;
    end
end
end
