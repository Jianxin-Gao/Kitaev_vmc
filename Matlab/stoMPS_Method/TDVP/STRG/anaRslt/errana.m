load('1.mat')
Bd = 5;
% semilogy(anaHn.err(Bd,:), '*')
% xlabel('n', 'FontSize', 20)
% ylabel('Trunction error', 'FontSize', 20)
% title({'H^n |CPS\rangle', ['\chi = ', num2str(Para.MCrit), '; L = ', num2str(Para.L)]}, 'FontSize', 20)
% keyboard

plot(anaHn.entropy(Bd,:), '*')
xlabel('n', 'FontSize', 20)
ylabel('Entanglement entropy', 'FontSize', 20)
title({'H^n |CPS\rangle', ['\chi = ', num2str(Para.MCrit), '; L = ', num2str(Para.L), '; \beta = ', num2str(Para.beta)]}, 'FontSize', 20)
keyboard


% semilogy(anaM.err(Bd,:), '*')
% xlabel('n', 'FontSize', 20)
% ylabel('Trunction error', 'FontSize', 20)
% title({'|METTS\rangle with n order expansion', ['D = ', num2str(Para.MCrit), '; L = ', num2str(Para.L)]}, 'FontSize', 20)
% keyboard

plot(anaM.entropy(Bd,:), '*')
xlabel('n', 'FontSize', 20)
ylabel('Entanglement entropy', 'FontSize', 20)
title({'|METTS\rangle with n order expansion', ['\chi = ', num2str(Para.MCrit), '; L = ', num2str(Para.L), '; \beta = ', num2str(Para.beta)]}, 'FontSize', 20)
keyboard
