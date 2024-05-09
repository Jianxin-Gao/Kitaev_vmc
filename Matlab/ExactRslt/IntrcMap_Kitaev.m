function [ Intr ] = IntrcMap_Kitaev(Model_Para)

L = Model_Para.L;

int_num = (L - 1) * L * 3 / 2;
int_cell = cell(int_num, 4);
count = 1;

% ========for Sx interaction==========
for col = 0:1:L-2
    if mod(col, 2) == 0
        for row = 0:2:L-2
            int_cell{count, 1} = L*col+row+1;
            int_cell{count, 2} = L*col+row+1 + L;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = Model_Para.Jx;
            count = count + 1;
        end
    else
        for row = 1:2:L-1
            int_cell{count, 1} = L*col+row+1;
            int_cell{count, 2} = L*col+row+1 + L;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = Model_Para.Jx;
            count = count + 1;
        end
    end
end
% ========for Sy interaction==========
for col = 0:1:L-1
    if mod(col, 2) == 0
        for row = 1:2:L-2
            int_cell{count, 1} = L*col+row+1 + 1;
            int_cell{count, 2} = L*col+row+1;
            int_cell{count, 3} = 'Sy';
            int_cell{count, 4} = Model_Para.Jy;
            count = count + 1;
        end
    else
        for row = 0:2:L-2
            int_cell{count, 1} = L*col+row+1 + 1;
            int_cell{count, 2} = L*col+row+1;
            int_cell{count, 3} = 'Sy';
            int_cell{count, 4} = Model_Para.Jy;
            count = count + 1;
        end
    end
end
% ========for Sz interaction==========
for col = 0:1:L-1
    if mod(col, 2) == 0
        for row = 0:2:L-2
            int_cell{count, 1} = L*col+row+1;
            int_cell{count, 2} = L*col+row+1 + 1;
            int_cell{count, 3} = 'Sz';
            int_cell{count, 4} = Model_Para.Jz;
            count = count + 1;
        end
    else
        for row = 1:2:L-2
            int_cell{count, 1} = L*col+row+1;
            int_cell{count, 2} = L*col+row+1 + 1;
            int_cell{count, 3} = 'Sz';
            int_cell{count, 4} = Model_Para.Jz;
            count = count + 1;
        end
    end
end

Intr = struct('Index1', int_cell(:,1), 'Index2', int_cell(:,2), ...
    'Intr_Type', int_cell(:,3),  'CS', int_cell(:,4));

