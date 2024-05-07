HnE = zeros(1, 45);
ME = zeros(1, 45);
Bd = 5;
for i = 1:1:100
    load([num2str(i), '.mat']);
    HnE = HnE + anaHn.entropy(Bd,1:1:45);
    ME = ME + anaM.entropy(Bd,1:1:45);
end
HnE = HnE/100;
ME = ME/100;

D1 = Para.MCrit;

load('SETTNRES')
D2 = Para.MCrit;

plot(HnE, '*'); hold on
plot(ME, 'o'); hold on
plot(anaHn.entropy(Bd,:), '*'); hold on
plot(anaM.entropy(Bd,:), 'o'); hold on
set(gca, 'FontSize', 15);
set(gca, 'LineWidth', 1.5)
xlabel('n', 'FontSize', 20)
ylabel('Entanglement entropy', 'FontSize', 20)
title({['Sample number = 100', '; L = ', num2str(Para.L), '; \beta = ', num2str(Para.beta)]}, 'FontSize', 20)
legend({'CompMPS', 'METTS', 'CompMPO', 'SETTN'}, 'FontSize', 20, 'Location', 'northwest')
axis([0 45 0 1.5])
% saveas(gcf, 'fig/EE.png')
hold off 