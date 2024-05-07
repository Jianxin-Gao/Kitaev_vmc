function [ H ] = MPS2V( MPO )

len = length(MPO);
T = MPO{1};
for i = 2:1:len
    T = contract(MPO{i}, 1, T, 1);
end
% T = permute(T, [1:2:(2 * len - 1), 2:2:(2 * len)]);
H = reshape(T, [2^len, 1]);
end

