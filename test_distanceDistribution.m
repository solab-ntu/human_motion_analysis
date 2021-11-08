
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

id1 = 6;
id2 = 8;
mu1 = load_template(template_names(id1));
mu2 = load_template(template_names(id2));

D1 = struct();
D2 = struct();
for i = 1:10
    load("Templates\\Raw_" + template_names(i) + ".mat")
    
    dist1 = zeros(1, length(IK));
    dist2 = zeros(1, length(IK));
    for j = 1:length(IK)
        data = IK(j).data;
        d1 = [];
        d2 = [];
        
        if size(data, 1) > size(mu1, 1)
            for k = 0:(size(data, 1) - size(mu1, 1))
                d1 = [d1, myMSM(data((1:size(mu1, 1))+k, :), mu1)];
            end
        elseif size(data, 1) < size(mu1, 1)
            for k = 0:(size(mu1, 1) - size(data, 1))
                d1 = [d1, myMSM(data, mu1((1:size(data, 1))+k, :))];
            end
        else
            d1 = myMSM(data, mu1);
        end

        if size(data, 1) > size(mu2, 1)
            for k = 0:(size(data, 1) - size(mu2, 1))
                d2 = [d2, myMSM(data((1:size(mu2, 1))+k, :), mu2)];
            end
        elseif size(data, 1) < size(mu2, 1)
            for k = 0:(size(mu2, 1) - size(data, 1))
                d2 = [d2, myMSM(data, mu2((1:size(data, 1))+k, :))];
            end
        else
            d2 = myMSM(data, mu2);
        end
        
        dist1(j) = mean(d1);
        dist2(j) = mean(d2);
    end
    
    D1(i).dist = dist1;
    D2(i).dist = dist2;
end

%%

names = ["Stand", "Squat Hd", "Squat Dn", "Squat Up", "Walk FL", ...
         "Walk FR", "Stair AL", "Stair AR", "Stair DL", "Stair DR"];

figure()
hold on
for i = 1:10
    plot(D1(i).dist, D2(i).dist, '.', 'markerSize', 8)
end
hold off
xlabel(join(split(template_names(id1),"_"), 1))
ylabel(join(split(template_names(id2),"_"), 1))
legend(names)
box on


