function [ tau_step, beta_c ] = GetStepFunc(tau, cric, max_beta)
Step = floor(log(cric(1)/tau)/log(2));
beta_c = tau * 2.^(2:Step);
tau_step = tau * 2.^(1:Step-1)/4;
beta_c = beta_c * 0.999;
beta_c = [beta_c, cric(1), max_beta];
tau_step = [tau_step, (cric(1)-tau*2.^Step)/4, cric(2:end)];
end