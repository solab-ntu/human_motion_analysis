
template_names = ["Stand";
                  "Squat_Hold";
                  "Squat_Down";
                  "Squat_Up";
                  "Walk_FL";
                  "Walk_FR";
                  "Stair_AL";
                  "Stair_AR";
                  "Stair_DL";
                  "Stair_DR"];

id = 10;
mu = load_template(template_names(id));
temp = load_mfpca(template_names);

D1 = struct();
D2 = struct();
D3 = struct();
D4 = struct();
D5 = struct();
for i = 1:size(template_names, 1)
    load("Templates\\Raw_" + template_names(i) + ".mat")
    
    dist1 = zeros(1, length(IK));
    dist2 = zeros(1, length(IK));
    dist3 = zeros(1, length(IK));
    dist4 = zeros(1, length(IK));
    dist5 = zeros(1, length(IK));
    for j = 1:length(IK)
        data = IK(j).data;
        d1 = [];
        d2 = [];
        d3 = [];
        d4 = [];
        if size(data, 1) > size(mu, 1)
            for k = 0:(size(data, 1) - size(mu, 1))
                d1 = [d1, myED(data((1:size(mu, 1))+k, :), mu)];
                d2 = [d2, myDTW(data((1:size(mu, 1))+k, :), mu)];
                d3 = [d3, myMSM(data((1:size(mu, 1))+k, :), mu)];
                d4 = [d4, myElasticDist(data((1:size(mu, 1))+k, :), mu)];
            end
        elseif size(data, 1) < size(mu, 1)
            for k = 0:(size(mu, 1) - size(data, 1))
                d1 = [d1, myED(data, mu((1:size(data, 1))+k, :))];
                d2 = [d2, myDTW(data, mu((1:size(data, 1))+k, :))];
                d3 = [d3, myMSM(data, mu((1:size(data, 1))+k, :))];
                d4 = [d4, myElasticDist(data, mu((1:size(data, 1))+k, :))];
            end
        else
            d1 = myED(data, mu);
            d2 = myDTW(data, mu);
            d3 = myMSM(data, mu);
            d4 = myElasticDist(data, mu);
        end

        dist1(j) = mean(d1);
        dist2(j) = mean(d2);
        dist3(j) = mean(d3);
        dist4(j) = mean(d4);
        dist5(j) = mfpcaDist(data, temp(id));
    end
    
    D1(i).dist = sort(dist1);
    D2(i).dist = sort(dist2);
    D3(i).dist = sort(dist3);
    D4(i).dist = sort(dist4);
    D5(i).dist = sort(dist5);
    
    if i == id
        disp('   ED        DTW       MSM       FRM       mFPCA')
        disp([sqrt(mean(dist1.^2)), sqrt(mean(dist2.^2)), sqrt(mean(dist3.^2)), sqrt(mean(dist4.^2)), sqrt(mean(dist5.^2))])
    end
end


ranges = [250, 220, 150, 150, 200, 200, 180, 180, 200, 200];
ranges3 = [80, 60, 40, 40, 60, 60, 60, 60, 60, 60];
ranges5 = [50, 20, 40, 40, 60, 60, 80, 80, 60, 60];

figure()
for i = 1:10
    w = 2;
    
    subplot(10, 5, 5*i-4)
    hold on
    histogram(D1(i).dist, 'BinWidth' , w, 'BinLimits', [0 ranges(id)], 'LineWidth', 0.1);
    yticklabels({})
    hold off
    grid on
    grid minor
    box on
    if i < 10
        xticklabels({})
    end

    subplot(10, 5, 5*i-3)
    hold on
    histogram(D2(i).dist, 'BinWidth' , w, 'BinLimits', [0 ranges(id)], 'LineWidth', 0.1);
    yticklabels({})
    hold off
    grid on
    grid minor
    box on
    if i < 10
        xticklabels({})
    end
    
    subplot(10, 5, 5*i-2)
    hold on
    histogram(D3(i).dist, 'BinWidth' , 1, 'BinLimits', [0 ranges3(id)], 'LineWidth', 0.1);
    yticklabels({})
    hold off
    grid on
    grid minor
    box on
    if i < 10
        xticklabels({})
    end
    
    subplot(10, 5, 5*i-1)
    hold on
    histogram(D4(i).dist, 'BinWidth' , w, 'BinLimits', [0 ranges(id)], 'LineWidth', 0.1);
    yticklabels({})
    hold off
    grid on
    grid minor
    box on
    if i < 10
        xticklabels({})
    end
    
    subplot(10, 5, 5*i)
    hold on
    histogram(D5(i).dist, 'BinWidth' , 1, 'BinLimits', [0 ranges5(id)], 'LineWidth', 0.1);
    yticklabels({})
    hold off
    grid on
    grid minor
    box on
    if i < 10
        xticklabels({})
    end
end

%%

function d = mfpcaDist(X, temp)
    
    t = (0:1/(size(X, 1) - 1):1)';
    tq = (0:1/(size(temp.mu, 1) - 1):1)';
    X = interp1(t, X, tq, 'makima');
    
    X = (X - temp.mu)./temp.std;
    k = length(temp.lambda);
    score = zeros(1, k);
    
    for i = 1:k
        score(i) = sum(sum(X.*temp.eig_func(:, :, i)));
    end
    d = sqrt(sum(score'.^2./temp.lambda));
end

function temp = load_mfpca(template_names)
    
    temp = struct();
    
    for i = 1:size(template_names, 1)
        
        fin = load("Templates\\Aligned_" + template_names(i) + ".mat");
        n = length(fin.IK); % number of samples
        l = size(fin.IK(1).data, 1);
        
        Data = zeros(l, n, 6);
        mu = zeros(l, 6);
        std = zeros(l, 6);
        for j = 1:n
            Data(:, j, 1) = fin.IK(j).data(:, 1);
            Data(:, j, 2) = fin.IK(j).data(:, 2);
            Data(:, j, 3) = fin.IK(j).data(:, 3);
            Data(:, j, 4) = fin.IK(j).data(:, 4);
            Data(:, j, 5) = fin.IK(j).data(:, 5);
            Data(:, j, 6) = fin.IK(j).data(:, 6);
        end
        mu(:, 1) = mean(Data(:, j, 1), 2);
        mu(:, 2) = mean(Data(:, j, 2), 2);
        mu(:, 3) = mean(Data(:, j, 3), 2);
        mu(:, 4) = mean(Data(:, j, 4), 2);
        mu(:, 5) = mean(Data(:, j, 5), 2);
        mu(:, 6) = mean(Data(:, j, 6), 2);
        std(:, 1) = sqrt(mean((Data(:, :, 1) - repmat(mu(:, 1), 1, n)).^2, 2));
        std(:, 2) = sqrt(mean((Data(:, :, 2) - repmat(mu(:, 2), 1, n)).^2, 2));
        std(:, 3) = sqrt(mean((Data(:, :, 3) - repmat(mu(:, 3), 1, n)).^2, 2));
        std(:, 4) = sqrt(mean((Data(:, :, 4) - repmat(mu(:, 4), 1, n)).^2, 2));
        std(:, 5) = sqrt(mean((Data(:, :, 5) - repmat(mu(:, 5), 1, n)).^2, 2));
        std(:, 6) = sqrt(mean((Data(:, :, 6) - repmat(mu(:, 6), 1, n)).^2, 2));
        
        C = [];
        for j = 1:6
           Cj = [];
           for k = 1:6
               Cjk = ((Data(:, :, j) - repmat(mu(:, j), 1, n))./repmat(std(:, j), 1, n))*((Data(:, :, k) - repmat(mu(:, k), 1, n))./repmat(std(:, k), 1, n))';
               Cj = [Cj, Cjk/n];
           end
           C = [C; Cj];
        end

        [U, V, ~] = svd(C);
        V = diag(V);
        k = 1;
        while sum(V(1:k)) < sum(V)*0.9
            k = k + 1;
        end

        eig_func = zeros(l, 6, k);
        for j = 1:k
            eig_func(:, :, j) = reshape(U(:,j), [l, 6]);
        end
        
        temp(i).mu = mu;
        temp(i).std = std;
        temp(i).eig_func = eig_func;
        temp(i).lambda = V(1:k);
    end
end

