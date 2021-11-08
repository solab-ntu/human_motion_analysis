function [mean, std, time] = load_template(fname)
    
    if ~isfile("Templates\\" + fname + ".mat")
        disp("Can't find " + fname);
        return
    end
    
    fin = load("Templates\\" + fname + ".mat");
    % n = numel(fieldnames(fin)); n must be 13
    
    time = fin.time;
    mean = zeros(length(time), 6);
    std = zeros(length(time), 6);
    
    mean(:, 1) = fin.hip_angle_l;
    mean(:, 2) = fin.hip_angle_r;
    mean(:, 3) = fin.knee_angle_l;
    mean(:, 4) = fin.knee_angle_r;
    mean(:, 5) = fin.ankle_angle_l;
    mean(:, 6) = fin.ankle_angle_r;
    
    std(:, 1) = fin.hip_angle_l_std;
    std(:, 2) = fin.hip_angle_r_std;
    std(:, 3) = fin.knee_angle_l_std;
    std(:, 4) = fin.knee_angle_r_std;
    std(:, 5) = fin.ankle_angle_l_std;
    std(:, 6) = fin.ankle_angle_r_std;
end