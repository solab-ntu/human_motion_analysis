function [d, Y_re, wp] = myMSM(Y, X)
% Move-Split-Merge refer to doi: 10.1109/TKDE.2012.88.
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
        
        if constraint_violate(i, j, 3)
            % Itakura parallelogram
            C(i, j) = Inf;
            continue
        end
        
        if i == 1 && j == 1
            C(1, 1) = dist(Y(1,:), X(1,:));
        elseif i == 1
            C(1, j) = C(1, j-1) + cost(X(1,:), Y(j,:), Y(j-1,:));
        elseif j == 1
            C(i, 1) = C(i-1, 1) + cost(Y(1,:), X(i,:), X(i-1,:));
        else
            d = [C(i-1, j-1) + dist(Y(j,:), X(i,:)),...
                 C(i, j-1) + cost(X(i,:), Y(j,:), Y(j-1,:)),...
                 C(i-1, j) + cost(Y(j,:), X(i,:), X(i-1,:))];
            C(i, j) = min(d);
        end
    end
end
d = sqrt(C(nx, ny)/min([nx, ny]));
% d = C(nx, ny)/min([nx, ny]);

% if number of output is 2, then find aligned Y_re
if nargout >= 2
    [ix, iy] = traceback(C);

    P = [0, 0];
    for i = 1:(length(ix)-1)
        if ix(i) == ix(i+1)
            continue
        else
            P = [P; ix(i) iy(i)];
        end
    end
    P = [P; ix(end) iy(end)];
    P = [P(:,1)/nx, P(:,2)/ny];
    
    wp = makima(P(:,1), P(:,2), linspace(0,1,nx));
    Y_re = interp1(linspace(0,1,ny)', Y, wp', 'makima');
end

% ================= functions ====================
function d = dist(Y, X) % metric
% one row of data
    d = sum((Y - X).^2);
    % d = sum(abs(Y - X));
end 

function c = cost(Y1, X1, X0)
% c0 is a system parameter for cost of split and merge
    c0 = 1;
    if nargin == 5
        ds = [dist(Y1, X1), dist(X1, X0)];
        if max(ds) <=  dist(Y1, X0)
            c = c0;
        else
            c = c0 + min(ds);
        end
    elseif nargin == 4
        ds = [dist(Y1, X1), dist(X1, X0)];
        if max(ds) <=  dist(Y1, X0)
            c = c0;
        else
            c = c0 + min(ds);
        end
    end
end

function frag = constraint_violate(i, j, slope)
% parallelogram
    if nargin == 2
        slope = 3;
    end
    
    if slope == 0
        frag = false;
        return
    end
    
    Smax = max([slope*ny/nx, 1]);
    Smin = min([1/slope*ny/nx, 1]);
    
    f = @(x, y, p1, p2) (Smax*(x-p1)-(y-p2))*(Smin*(x-p1)-(y-p2));
    
    if f(i, j, 1, 1) <= 0 && f(i, j, nx, ny) <= 0
        frag = false;
    else
        frag = true;
    end
end


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