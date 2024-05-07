HnE = zeros(100, 45);
ME = zeros(100, 45);
Bd = 5;
count = 0;
for i = 1:1:100
    load([num2str(i), '.mat']);
    count = count + 1;
    HnE(count,:) = anaHn.entropy(Bd,1:1:45);
    ME(count,:) = anaM.entropy(Bd,1:1:45);
    % plot(1:1:45, anaM.entropy(Bd,1:1:45), '*'); hold on
end

HnE = HnE(1:1:count, :);
ME = ME(1:1:count, :);
sigmaHnE = std(HnE, 1);
meanHnE = mean(HnE);

sigmaME = std(ME, 1);
meanME = mean(ME);

errorbar(meanHnE,sigmaHnE, 'o', 'MarkerSize',10, 'MarkerFaceColor', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on 
errorbar(meanME,sigmaME, 'o', 'MarkerSize',10, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)

load('SETTNRES')
D2 = Para.MCrit;

plot(anaHn.entropy(Bd,:), '*', 'MarkerSize',10, 'LineWidth', 1.5); hold on
plot(anaM.entropy(Bd,:), '*', 'MarkerSize',10, 'LineWidth', 1.5); hold on
set(gca, 'FontSize', 15);
set(gca, 'LineWidth', 1.5)
xlabel('n', 'FontSize', 20)
ylabel('Entanglement entropy', 'FontSize', 20)
title({['Sample number = 100', '; L = ', num2str(Para.L), '; \beta = ', num2str(Para.beta)]}, 'FontSize', 20)
legend({'CompMPS', 'METTS', 'CompMPO', 'SETTN'}, 'FontSize', 20, 'Location', 'northwest')
axis([0 45 0 1.5])
saveas(gcf, 'fig/EE.png')
hold off 

% plot(sigmaME)