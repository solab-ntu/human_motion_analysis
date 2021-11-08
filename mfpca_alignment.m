
% used with mfpca_raw to plot

% load("Templates\\Raw_Squat_Down.mat");
% load("Templates\\Raw_Walk_FL.mat");
% load("Templates\\Raw_Stair_AL.mat");
load("Templates\\Raw_Stair_DL.mat");

% time warping alignment
% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

min_dist = Inf;
for i = 1:length(IK)
    d = 0;
    for j = 1:length(IK)
        if i ~= j
            % d = d + myED(IK(j).data, IK(i).data);
            % d = d + myDTW(IK(j).data, IK(i).data);
            % d = d + myMSM(IK(j).data, IK(i).data);
            d = d + myElasticDist(IK(j).data, IK(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_avg = IK(min_id).data;
id = min_id;

%% 

pre = Inf;
while true
    IK_align = struct();
    tmp = 0;
    
    for i = 1:length(IK)
        if i == id
            IK_align(i).data = IK_avg;
            id = 0;
        else
            % [d, data, wp] = myED(IK(i).data, IK_avg);
            % [d, data, wp] = myDTW(IK(i).data, IK_avg);
            [d, data, wp] = myMSM(IK(i).data, IK_avg);
            IK_align(i).data = data;
            IK_align(i).wpf = wp;
            tmp = tmp + d;
        end
    end

    IK_avg = IK_align(1).data;
    for i = 2:length(IK)
        IK_avg = IK_avg + IK_align(i).data;
    end
    IK_avg = IK_avg/length(IK);
    
    if abs(tmp - pre)/tmp < 0.01
        break
    else
        pre = tmp;
    end
end

IK = IK_align;

%% 

l = size(IK_avg, 1);

% initialize input curve
IK_beta = zeros(size(IK(1).data, 2), l, length(IK));
for i = 1:length(IK)
    Data = IK(i).data;
    Data = interp1(linspace(0,1,size(IK(i).data, 1))', Data, linspace(0,1,l)', 'makima');
    IK_beta(:,:,i) = Data';
end

% computing aligned curves
ESA_obj = fdacurve(IK_beta, false, l);
ESA_objn = ESA_obj.karcher_mean();


IK_ESA = struct();
for i = 1:length(IK)
    IK_ESA(i).wpf = ESA_objn.gams(:,i);
    IK_ESA(i).data = interp1(linspace(0,1,size(IK_beta(:,:,i),2)),IK_beta(:,:,i)',IK_ESA(i).wpf,'makima');
end

IK = IK_ESA;

%%

for i = 1:length(IK)
    data = [IK(i).data(:,2), IK(i).data(:,1), IK(i).data(:,4), IK(i).data(:,3), IK(i).data(:,6), IK(i).data(:,5)];
    IK(i).data = data;
end

%%
for i = 1:length(IK)
    IK(i).data = flip(IK(i).data, 1);
end

