function [] = plot_similarity(templates, duration, title_string) % ====
    
    template_names = {};
    for i = 1:length(templates)
        template_names = [template_names, templates(i).name];
    end
    
    l = length(templates(1).time);
    dl = max([round(l/80), 1]);
    mk = ['o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];
    mks = mk;
    while length(mks) < length(templates)
        mks = [mks, mk];
    end
    
    figure()
    box on
    hold on
    for i = 1:length(templates)
        plot(templates(i).time, templates(i).smlr, 'LineWidth', 1, 'Marker', mks(i), 'MarkerIndices', 1:dl:l, 'MarkerSize', 5)
    end
    hold off
    legend(template_names)
    title(title_string)
    xlim(duration)
    % ylim([0 1])
    xlabel("Time(sec)")
    grid on
    grid minor
end