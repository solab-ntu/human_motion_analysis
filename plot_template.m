function [] = plot_template(fname, title_str, yl)
% fname  file name
% tp     file type

    if nargin == 1
        title_str = join(split(fname,"_"), 1);
    end

    [M, S, t] = load_template(fname);
    
    % plot
    figure()
    subplot(1,2,1)
    box on
    hold on
    fill([t', fliplr(t')], [M(:,1)' + S(:,1)', fliplr(M(:,1)' - S(:,1)')], 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    fill([t', fliplr(t')], [M(:,3)' + S(:,3)', fliplr(M(:,3)' - S(:,3)')], 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    fill([t', fliplr(t')], [M(:,5)' + S(:,5)', fliplr(M(:,5)' - S(:,5)')], 1, 'FaceColor', '#EDB120', 'FaceAlpha', 0.2, 'LineStyle',' none');
    h1 = plot(t, M(:,1), 'Color', '#0072BD', 'LineWidth', 1);
    h2 = plot(t, M(:,3), 'Color', '#D95319', 'LineWidth', 1);
    h3 = plot(t, M(:,5), 'Color', '#EDB120', 'LineWidth', 1);
    hold off
    legend([h1, h2, h3], {'Hip', 'Knee', 'Ankle'})
    title(title_str + " (Left Leg)")
    xlabel("Time(sec)")
    xlim([t(1), t(end)])
    ylabel("Angle(degree)")
    grid on
    grid minor
    if nargin == 3
        ylim(yl)
    end
    
    subplot(1,2,2)
    box on
    hold on
    fill([t', fliplr(t')], [M(:,2)' + S(:,2)', fliplr(M(:,2)' - S(:,2)')], 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    fill([t', fliplr(t')], [M(:,4)' + S(:,4)', fliplr(M(:,4)' - S(:,4)')], 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    fill([t', fliplr(t')], [M(:,6)' + S(:,6)', fliplr(M(:,6)' - S(:,6)')], 1, 'FaceColor', '#EDB120', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    h1 = plot(t, M(:,2), 'Color', '#0072BD', 'LineWidth', 1);
    h2 = plot(t, M(:,4), 'Color', '#D95319', 'LineWidth', 1);
    h3 = plot(t, M(:,6), 'Color', '#EDB120', 'LineWidth', 1);
    hold off
    legend([h1, h2, h3], {'Hip', 'Knee', 'Ankle'})
    title(title_str + " (Right Leg)")
    xlabel("Time(sec)")
    xlim([t(1), t(end)])
    ylabel("Angle(degree)")
    grid on
    grid minor
    if nargin == 3
        ylim(yl)
    end
end