load('1.mat')
Bd = 5;

plot(anaHn.entropy(Bd,:), '*'); hold on
xlabel('n', 'FontSize', 20)
ylabel('Entanglement entropy', 'FontSize', 20)
title({['\chi = ', num2str(Para.MCrit), '; L = ', num2str(Para.L), '; \beta = ', num2str(Para.beta)]}, 'FontSize', 20)

plot(anaM.entropy(Bd,:), '*')
xlabel('n', 'FontSize', 20)
ylabel('Entanglement entropy', 'FontSize', 20)

legend({'CompMPS', 'METTS'}, 'FontSize', 20, 'Location', 'northwest')
hold off