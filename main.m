
template_names = ["Stand", "Stand";
          "Squat Dn", "Squat_Down";
          "Squat Hd", "Squat_Hold";
          "Squat Up",   "Squat_Up";
          "Walk FL",     "Walk_FL";
          "Walk FR",     "Walk_FR";
          "Stair AL",   "Stair_AL";
          "Stair AR",   "Stair_AR";
          "Stair DL",   "Stair_DL";
          "Stair DR",   "Stair_DR"];

fin = "Scenario_1";
load("Scenario\\" + fin + ".mat");
templates = initial_temp(template_names, size(scenario, 1));

%%

metric = 'ED';
templates = similarity(scenario, templates, metric);
plot_event(templates, metric)
% plot_distance(templates, 1)
% save("Result\\" + fin + "_" + metric + ".mat", 'templates')
% disp(fin + "_" + metric)

%%

metric = 'DTW';
fin = "Scenario_1_wgn3";
load("Result\\" + fin + "_" + metric + ".mat");
plot_event(templates, metric)
disp(fin + "  " +metric)
% title(join(split(fin, "_"), 1) + " " + metric)

%%

for fin = ["Scenario_1", "Scenario_1_Raw", "Scenario_1_wgn3", "Scenario_2"]
    % load("Scenario\\" + fin + ".mat");
    for metric = ["ED", "DTW", "MSM", "FRM", "PCA"]
        %templates = similarity(scenario, templates, metric);
        %save("Result\\" + fin + "_" + metric + ".mat", 'templates')
        %disp(fin + "_" + metric)
        
        load("Result\\" + fin + "_" + metric + ".mat");
        plot_event(templates, metric)
        if metric == "PCA"
            saveas(gcf, "figure\\6-3_" + fin + "_mFPCA.png")
        else
            saveas(gcf, "figure\\6-3_" + fin + "_" + metric + ".png")
        end
        disp(fin + "  " +metric)
        
    end
end

%%

function templates = initial_temp(template_names, len)
    
    templates = struct();
    
    for i = 1:size(template_names, 1)
        
        fin = load("Templates\\Raw_" + template_names(i, 2) + ".mat");
        % fin = load("Templates\\Aligned_" + template_names(i, 2) + ".mat");
        n = length(fin.IK); % number of samples
        
        l = 0;
        for j = 1:n
            l = l + size(fin.IK(j).data, 1);
        end
        l = round(l/n); % mean length of samples
        
        Data = zeros(l, n, 6);
        mu = zeros(l, 6);
        std = zeros(l, 6);
        for j = 1:n
            DataIn = fin.IK(j).data;
            t = 0:1/(size(DataIn, 1) - 1):1;
            tq = 0:1/(l - 1):1;
            DataIn = interp1(t', DataIn, tq', 'makima');
            Data(:, j, 1) = DataIn(:, 1);
            Data(:, j, 2) = DataIn(:, 2);
            Data(:, j, 3) = DataIn(:, 3);
            Data(:, j, 4) = DataIn(:, 4);
            Data(:, j, 5) = DataIn(:, 5);
            Data(:, j, 6) = DataIn(:, 6);
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
        
        templates(i).name = template_names(i, 1);
        templates(i).mean = mu;
        templates(i).std = std;
        templates(i).eig_func = eig_func;
        templates(i).lambda = V(1:k);
        templates(i).dist = zeros(len, 1);
    end
end