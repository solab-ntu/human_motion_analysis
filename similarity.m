function templates = similarity(scenario, templates, metric)

% switch metric
%     case 'ED'
%         sigma = [8.7027, 19.3944, 22.8257, 22.8257, 14.1550, 14.1550, 14.8650, 14.8650, 19.1579, 19.1579];
%     case 'DTW'
%         sigma = [8.7027, 19.3944, 9.8670, 9.8670, 8.5889, 8.5889, 10.3837, 10.3837, 11.3342, 11.3342];
%     case 'MSM'
%         sigma = [2.7578, 5.5523, 2.7240, 2.9103, 4.1388, 4.1388, 4.0975, 4.0975, 5.0983, 5.0983];
%     case 'FRM'
%         sigma = [8.6669, 19.6077, 19.2173, 18.5376, 9.5250, 9.5250, 11.7324, 11.7324, 12.4458, 12.4458];
%     case 'PCA'
%         sigma = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
% end

for t = 0:(size(scenario, 1)-1) % move window
   
    for i = 1:length(templates)
        
        dv = size(templates(i).mean, 1);
        tt = t + ((1-floor(dv/2)):ceil(dv/2))';
        tt = tt(tt >= 1 & tt <= length(scenario(:,1)));
        
        y_samp = scenario(tt, :);
        
        switch metric
            case 'ED'
                dist = myED(y_samp, templates(i).mean);
            case 'DTW'
                dist = myDTW(y_samp, templates(i).mean);
            case 'MSM'
                dist = myMSM(y_samp, templates(i).mean);
            case 'FRM'
                dist = myElasticDist(y_samp, templates(i).mean);
            case 'PCA'
                dist = mfpcaDist(y_samp, templates(i));
            case 'DTW+PCA'
                dist = mfpcaDTW(y_samp, templates(i));
        end
        
        templates(i).dist(t+1) = dist;
    end
    
%     if mod(t,100) == 0
%         disp(t)
%     end
end

function d = mfpcaDist(samp, temp)
    
    x = (0:1/(size(samp, 1) - 1):1)';
    xq = (0:1/(size(temp.mean, 1) - 1):1)';
    samp = interp1(x, samp, xq, 'makima');
    
    samp = (samp - temp.mean)./temp.std;
    k = length(temp.lambda);
    score = zeros(1, k);
    
    for j = 1:k
        score(j) = sum(sum(samp.*temp.eig_func(:, :, j)));
    end
    d = sqrt(sum(score'.^2./temp.lambda));
end

function d = mfpcaDTW(samp, temp)
    
    x = (0:1/(size(samp, 1) - 1):1)';
    xq = (0:1/(size(temp.mean, 1) - 1):1)';
    samp = interp1(x, samp, xq, 'makima');
    
    [~, samp, ~] = myDTW(samp, temp.mean);
    
    samp = (samp - temp.mean)./temp.std;
    k = length(temp.lambda);
    score = zeros(1, k);
    
    for j = 1:k
        score(j) = sum(sum(samp.*temp.eig_func(:, :, j)));
    end
    d = sqrt(sum(score'.^2./temp.lambda));
end

end