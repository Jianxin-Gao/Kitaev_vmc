function [beta_list] = Get_beta_expTau(input1, varargin)
% input1 = [tau0, betamax1]
% varargin{1} = [linear tau, betamax2]

% Get beta_list = [4tau0, 8tau0, 16tau0, ..., betamax1]
tau0 = input1(1);
betamax1 = input1(2);
stepMax = floor(log2(betamax1/tau0));
beta_list = 2.^(2:stepMax)*tau0;
if ~isempty(varargin)
    % Get beta_list2 = [betamax1, betamax1+4tau1, ..., betamax2]
    temp = varargin{1};
    tau1 = temp(1);
    betamax2 = temp(2);
    beta_list2 = (beta_list(end)+4*tau1):4*tau1:betamax2;
    beta_list = [beta_list, beta_list2];
end
end