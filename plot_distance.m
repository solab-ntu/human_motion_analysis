function [] = plot_distance(templates, type)

if nargin == 1
    type = 1;
end

if type == 1
    
    template_names = {};
    for i = 1:length(templates)
        template_names = [template_names, templates(i).name];
    end
    
    l = size(templates(1).dist, 1);
    dl = max([round(l/80), 1]);
    time = 0.01*((1:l)'-1);
    mk = ['o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];
    mks = mk;
    while length(mks) < length(templates)
        mks = [mks, mk];
    end
    
    figure()
    hold on
    for i = 1:length(templates)
        plot(time, templates(i).dist/sigma(i), 'LineWidth', 1, 'Marker', mks(i), 'MarkerIndices', 1:dl:l, 'MarkerSize', 5)
    end
    hold off
    legend(template_names)
    xlim([time(1), time(end)])
    xticks(0:2:time(end))
    xlabel("time(sec)")
    box on
    grid on
    % grid minor

elseif type == 2
    
    l = size(templates(1).dist, 1);
    time = 0.01*((1:l)'-1);
    
    figure()
    for i = 1:length(templates)
        subplot(10, 1, i)
        plot(time, templates(i).dist, 'LineWidth', 1)
        xlim([time(1), time(end)])
        xticks(0:2:time(end))
        box on
        grid on
        
        if i == 10
            xlabel("time(sec)")
        end
    end
end
    
end