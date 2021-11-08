function dy = myElasticDist(Y, X)

dt = 0.01;
Xq = sign(diff(X,1)).*sqrt(abs(diff(X,1))/dt);
Yq = sign(diff(Y,1)).*sqrt(abs(diff(Y,1))/dt);

wp = srvf_wp(Yq, Xq);
wp = makima(linspace(0,1,length(wp)), wp, linspace(0,1,size(X,1)));
Yn = interp1(linspace(0,1,size(Y,1))', Y, wp', 'makima');
dy = sqrt(sum(sum((X - Yn).^2))/min([size(X, 1), size(Y, 1)]));

% dy = srvf_dist(Yq, Xq);

% ============================================================
function wp = srvf_wp(Yq, Xq)

nx = size(Xq, 1);
ny = size(Yq, 1);
C = zeros(nx, ny);

kk = zeros(nx, ny);
for i = 1:nx
   for j = 1:ny
       if i == 1 && j == 1
           C(1,1) = sum((Yq(j,:) - Xq(i,:)).^2);
       elseif i == 1
           C(1,j) = Inf;
       else
           d = sum((repmat(Xq(i,:),j,1) - sqrt((j-1):(-1):0)'*Yq(j,:)).^2, 2) + C(i-1,1:j)';
           [C(i,j), kk(i,j)] = min(d);
       end
   end

end

[ix, iy] = traceback(kk);

P = [0, 0; ix, iy];
P = [P(:,1)/(nx), P(:,2)/(ny)];

xx = linspace(0,1,nx+1);
wp = makima(P(:,1), P(:,2), xx(2:end));

end

function d = srvf_dist(Yq, Xq) % =================================

nx = size(Xq, 1);
ny = size(Yq, 1);
C = zeros(nx, ny);

kk = zeros(nx, ny);
for i = 1:nx
   for j = 1:ny
       if i == 1 && j == 1
           C(1,1) = sum((Yq(j,:) - Xq(i,:)).^2);
       elseif i == 1
           C(1,j) = Inf;
       else
           dmin = Inf;              
           k0 = j;
           for k = 1:j
               d = sum((Xq(i,:) - Yq(j,:)*sqrt(j-k)).^2) + C(i-1,k);
               if d < dmin
                   dmin = d;
                   k0 = k;
               end
           end
           C(i,j) = dmin;
           kk(i,j) = k0;
       end
   end

end

[ix, iy] = traceback(kk);

P = [0, 0; ix, iy];

xx = linspace(0,1,nx+1);
yy = linspace(0,1,ny+1);
wp = makima(P(:,1)/(nx), P(:,2)/(ny), xx(2:end));
Yqn = interp1(yy(2:end)', Yq, wp', 'makima');

d = sqrt(sum(sum((Xq - Yqn.*sqrt(diff(P(:,2), 1))).^2)));

end

function [ix_out,iy_out] = traceback(K) % =======================
    P = [];    
    [i, j] = size(K);
    while i > 0 
        P = [P; i, j];
        j = K(i,j);
        i = i - 1;
    end
    ix_out = flip(P(:,1),1);
    iy_out = flip(P(:,2),1);
end

end