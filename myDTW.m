function [d, Y_re, wp] = myDTW(Y, X)
% Dynamic Time Warping
% input X is fixed, and align Y to X
% S is stantard variation of X, asigned to 1 by default
% number of columns is numder of dimensions
% number of rows is numder of data point, namely the length of X, Y
% output Y_re is aligned Y, and d is dtw distance divided by mean length of X, Y

% example input
% X = [6 ,8 ,9 ,8 ,6 ,3 ,1 ,0 ,1, 3]';
% Y = [9 ,8 ,6 ,3 ,1 ,0 ,1, 3, 6, 8]';

nx = size(X,1);
ny = size(Y,1);

C = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        d = sum((X(i,:) - Y(j,:)).^2); % metric
        if i == 1 && j == 1
            C(1, 1) = d;
        elseif i == 1
            C(1, j) = d + C(1, j-1);
        elseif j == 1
            C(i, 1) = d + C(i-1, 1);
        else
            C(i, j) = d + min([C(i-1, j), C(i, j-1), C(i-1, j-1)]);
        end
    end
end
d = sqrt(C(nx, ny)/min([nx, ny]));

% using the matlab build-in dtw for testing
% [d,iy,ix] = dtw(Y',X');
% d = sqrt(d/min([size(X,1), size(Y,1)]));

% if number of output is 2, then find aligned Y_re

if nargout >= 2
    [ix, iy] = traceback(C);

%     Y_re = zeros(size(X)); 
%     s = 0;
%     count = 0;
%     for i = 1:length(ix)
%         s = s + Y(iy(i),:);
%         count = count + 1;
%         if i == length(iy)
%             Y_re(ix(i),:) = s/count;
%         elseif ix(i) ~= ix(i+1)
%             Y_re(ix(i),:) = s/count;
%             s = 0;
%             count = 0;
%         end
%     end
    
    P = [0, 0];
    for i = 1:(length(ix)-1)
        if ix(i) == ix(i+1)
            continue
        else
            P = [P; ix(i) iy(i)];
        end
    end
    P = [P; ix(end) iy(end)];
    
    xx = linspace(0,1,nx+1);
    yy = linspace(0,1,ny+1);
    wp = makima(P(:,1)/(nx), P(:,2)/(ny), xx(2:end));
    Y_re = interp1(yy(2:end)', Y, wp', 'makima');
end

% ================= functions ====================
function [ix_out,iy_out] = traceback(C)
    
    m = size(C,1);
    n = size(C,2);
    ix = zeros(m+n,1);
    iy = zeros(m+n,1);
    ix(1) = m;
    iy(1) = n;
    
    i = m;
    j = n;
    k = 1;
    while i>1 || j>1
        if j == 1
            i = i-1;
        elseif i == 1
            j = j-1;
        else
            cij = C(i-1,j-1);
            ci = C(i-1,j);
            cj = C(i,j-1);
            i = i - (ci<=cj | cij<=cj | cj~=cj);
            j = j - (cj<ci | cij<=ci | ci~=ci);
        end
        k = k+1;
        ix(k) = i;
        iy(k) = j;
    end

    ix_out = zeros(k,1);
    iy_out = zeros(k,1);
    for id = 1:k
        ix_out(id) = ix(k-id+1);
        iy_out(id) = iy(k-id+1);
    end
end
end


    