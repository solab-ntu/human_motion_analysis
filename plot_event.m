function [] = plot_event(templates, metric)
    
    w = 0.9; % width of bar
    l = size(templates(1).dist, 1);
    time = 0.01*((1:l)'-1);
    clr = [0 0.219 0.602]; % color

    template_names = {};
    for i = 1:length(templates)
        template_names = [template_names, templates(i).name];
    end
    len = length(template_names);
    
    switch metric
        case 'ED'
            sigma = [8.7027, 19.3944, 22.8257, 22.8257, 14.1550, 14.1550, 14.8650, 14.8650, 19.1579, 19.1579];
        case 'DTW'
            sigma = [8.7027, 19.3944, 9.8670, 9.8670, 8.5889, 8.5889, 10.3837, 10.3837, 11.3342, 11.3342];
        case 'MSM'
            sigma = [2.7578, 5.5523, 2.7240, 2.9103, 4.1388, 4.1388, 4.0975, 4.0975, 5.0983, 5.0983];
        case 'FRM'
            sigma = [8.6669, 19.6077, 19.2173, 18.5376, 9.5250, 9.5250, 11.7324, 11.7324, 12.4458, 12.4458];
        case 'PCA'
            sigma = [1.4142, 1.4142, 2.6708, 2.6708, 2.5449, 2.5449, 2.0506, 2.0506, 2.9992, 2.9992];
        otherwise
            sigma = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    end
    
    sftmax = [templates(1).dist, templates(2).dist, templates(3).dist, templates(4).dist, ...
              templates(5).dist, templates(6).dist, templates(7).dist, templates(8).dist, ...
              templates(9).dist, templates(10).dist];
    % sftmax = exp(- (sftmax./repmat(sigma, l, 1)).^2);
    % sftmax = 1./(1+(sftmax./repmat(3*sigma, l, 1)).^6);
    sftmax = exp(- sftmax)./repmat(sum(exp(- sftmax), 2), 1, 10)./(1+(sftmax./repmat(3*sigma, l, 1)).^6);
    
    % figure()
    clf
    
    box on
    hold on
    for i = 1:l
        for j = 1:10
            % p = [x, y, dx, dy] position and size of rectangle
            p = [time(i) - 0.005, len + 1 - j - w/2, 0.01, w];    
            rectangle('Position', p, 'FaceColor', [clr sftmax(i, j)], 'EdgeColor', 'none')
        end
    end
    xlim([time(1), time(end)])
    xticks(0:2:time(end))
    ylim([0, len+1])
    yticks(1:len)
    yticklabels(fliplr(template_names))
    xlabel("time(sec)")
    set(gca, 'xminorgrid', 'on', 'XMinorTick', 'on');
    grid on
    
    myMap = [linspace(1,clr(1),100)', linspace(1,clr(2),100)', linspace(1,clr(3),100)'];
    cb = colorbar;
    cb.Label.String = 'Score';
    colormap(myMap)
end