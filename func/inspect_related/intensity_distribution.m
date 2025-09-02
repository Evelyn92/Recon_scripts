figure;
for ind = 1:6
    subplot(3,2,ind)

    nb = 256;                         % # of bins (common choices: 64–512)
    xd = norm_single_volume(x{ind}(:));
    edges = linspace(min(xd), max(xd), nb+1);
    histogram(xd(:), edges, 'Normalization','pdf');     % probability density
    % xlim([0,0.5])
    xlabel('Intensity'); ylabel('Probability / bin');
    title(['Normalized intensity distribution: ', num2str(ind)]);
end