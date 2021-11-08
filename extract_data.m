% colunms of txt data contain:
% ["Time", "Pelvis", "Lumbar", "Hip R", "Hip L", "Knee R", "Knee L", "Ankle R", "Ankle L"]
% colunms of extracted data contain:
% ["Time", "Hip R", "Hip L", "Knee R", "Knee L", "Ankle R", "Ankle L"]

% note:
% colunms of new data of motion and scenario
% ["Hip L", "Hip R", "Knee L", "Knee R", "Ankle L", "Ankle R"]
% but i'm too lazy to correct this code:3

function [Data, DataformOut] = extract_data(fname, t1, t2, dt, save_output, s)
% read data from .ini and .txt which are coverted from trc, mot, sto
% data must be preprocessed by ExtractOpensimResult.py
% DataformOut = [nrow, ncol, isDegree, items...] in string
% Data cantains values of items each row


fin = fopen("Data_txt\\" + fname + ".ini", 'r');
Dataform = textscan(fin, "%s");
fclose(fin);

if nargout == 2
    DataformOut = [];
    for i = 1 : length(Dataform{1})
        DataformOut = [DataformOut string(Dataform{1}{i})];
    end
end

nRow = str2num(Dataform{1}{1});
nCol = str2num(Dataform{1}{2});

fin = fopen("Data_txt\\" + fname + ".txt", 'r');
format = "%f";
for i = 1 : (nCol - 1)
    format = format + "%f";
end
Data = fscanf(fin, format, [nCol nRow])';
fclose(fin);

% data extract and 1-D interpolation
% 'linear'(defuat), 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', 'spline'
% hip = hip + pelvis
% dt = 0.01;

Data = Data((Data(:,1) >= t1)&(Data(:,1) <= t2), :);

if size(Data,2) == 9
    Data(:,[4,5]) = Data(:,[4,5]) + Data(:,2);
    Data(:,[2,3]) = [];
elseif size(Data,2) == 8
    Data(:,[3,4]) = Data(:,[3,4]) + Data(:,2);
    Data(:,2) = [];
end

% Data(:,1) = Data(:,1) - Data(1,1);
tt = (0 :dt :(Data(end,1)-Data(1,1)))';
Data = interp1(Data(:,1) - Data(1,1), Data, tt, 'makima');

if nargin == 4 
    save_output = false;
end

if save_output
    if nargin == 5
        s = "";
    end
    time = Data(:,1);
    hip_angle_r = Data(:,2);
    hip_angle_l = Data(:,3);
    knee_angle_r = Data(:,4);
    knee_angle_l = Data(:,5);
    ankle_angle_r = Data(:,6);
    ankle_angle_l = Data(:,7);

    path = "Data_mat\\" + fname + s + ".mat";
    save(path, 'time', 'hip_angle_r', 'hip_angle_l', ...
               'knee_angle_r', 'knee_angle_l', ...
               'ankle_angle_r', 'ankle_angle_l');  
end

end