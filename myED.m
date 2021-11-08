function [dmin, Y_re, wp] = myED(Y, X)
% Euclidean Distance
% input X is fixed, and align Y to X by offset
% S is stantard variation of X, asigned to 1 by default
% number of columns is numder of dimensions
% number of rows is numder of data point, namely the length of X, Y
% output Y_re is aligned Y, and d is dtw distance divided by mean length of X, Y

% example input
% X = [6 ,8 ,9 ,8 ,6 ,3 ,1 ,0 ,1, 3]';
% Y = [9 ,8 ,6 ,3 ,1 ,0 ,1, 3, 6, 8]';

dmin = Inf;
offset = 0;
nx = size(X,1);
ny = size(Y,1);
m = size(X,2);

if nx >= ny % n(X) >= n(Y)
    for i = 1:(nx-ny+1)
        d = sum(sum((Y - X(i:(i+ny-1), :)).^2));
        if d < dmin
            dmin = d;
            offset = i - 1;
        end
    end
else         % n(X) < n(Y)
    for i = 1:(ny-nx+1)
        d = sum(sum((Y(i:(i+nx-1), :) - X).^2, 2));
        if d < dmin
            dmin = d;
            offset = i - 1;
        end
    end
end

dmin = sqrt(dmin/min([nx ny]));

% if number of output is 2, then find aligned Y_re
if nargout >= 2
    if nx >= ny
        ix = [1:(offset+ny), (ny+offset+1):nx];
        iy = [ones(1,offset), 1:ny, ones(1,nx-ny-offset)*ny];
    else
        ix = [ones(1,offset), 1:nx, ones(1,ny-nx-offset)*nx];
        iy = [1:(offset+nx), (nx+offset+1):ny];
    end
    
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

% % if number of output is 2, then find aligned Y_re
% if nargout >= 2
%     if nx >= ny
%         
%         if offset == 0
%             Y_re = Y;
%         else
%             Y_re = [ones(offset,1)*Y(1,:); Y];
%         end
%         
%         if length(Y_re') < nx
%             Y_re = [Y_re; ones(nx-size(Y_re,1),1)*Y(end,:)];
%         end
%     else
%         Y_re = Y;
%         
%         if offset ~= 0
%             Y_re(1:offset, :) = [];
%         end
%         
%         if size(Y_re,1) > nx
%             Y_re((nx+1):end, :)=[];
%         end
%     end
% end
% 
% % warping function
% if nargout == 3
%     if nx >= ny
%         if offset == 0
%             wp = 1:ny;
%         else
%             wp = [ones(1, offset), 1:ny];
%         end
%         
%         if length(wp') < nx
%             wp = [wp, ones(1, nx-length(wp))*ny];
%         end
%     else
%         wp = 1:ny;
%         
%         if offset ~= 0
%             wp(1:offset) = [];
%         end
%         
%         if length(wp) > nx
%             wp((nx+1):end)=[];
%         end
%         
%         wp(wp>nx) = nx;
%     end
% end
end