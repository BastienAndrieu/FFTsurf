clc; clear; close all


fid = fopen('../polylines/result.dat','r');
ninter = str2num(fgetl(fid));
ipls = zeros(ninter,2);
for i = 1:ninter
    ipls(i,:) = str2num(fgetl(fid));
end
lamb = zeros(ninter,2);
for i = 1:ninter
    lamb(i,:) = str2num(fgetl(fid));
end
fclose(fid);

for i = 1:2
    fid = fopen(sprintf('../polylines/xy%d.dat',i),'r');
    poly(i).np = str2num(fgetl(fid));
    poly(i).xy = zeros(poly(i).np,2);
    for j = 1:poly(i).np
        poly(i).xy(j,:) = str2num(fgetl(fid));
    end
    fclose(fid);
end


[~,iinter] = max(lamb(:,1));
ipls = ipls(iinter,:);
lamb = lamb(iinter,:);

uvtmp = poly(2).xy(ipls(2)+1,:);
duvtmp = poly(1).xy(2,:) - uvtmp;

tng = [
    -0.54243555491001216        4.7309090701704157
%     7.3188064138298969        1.4780358447249984 
    ];




figure;
hold on

a = poly(1).xy(1,:);
v = poly(1).xy(2,:) - a;

quiver(a(1), a(2), v(1), v(2), ...
    0, 'b', 'maxheadsize', 1);
for i = 2
    plot(poly(i).xy(:,1), poly(i).xy(:,2), '.-');
    plot(poly(i).xy(1,1), poly(i).xy(1,2), 'o');
    k = ipls(i);
    xy = (1.0 - lamb(i)) * poly(i).xy(k,:) + lamb(i) * poly(i).xy(k+1,:);
    plot(xy(1), xy(2), 'k*');
end

quiver(uvtmp(1), uvtmp(2), duvtmp(1), duvtmp(2), 0, 'r', 'maxheadsize', 1);


l = norm(v) / norm(tng);
quiver(uvtmp(1), uvtmp(2), l*tng(1), l*tng(2), 0, 'g', 'maxheadsize', 1);

axis image


p = poly(2).xy(ipls(2)+1,:) - poly(2).xy(ipls(2),:);
p = -p;
w = [-tng(2), tng(1)] / norm(tng);

r = 0.5 * sign(dot(tng,p)) * sum(p.^2) / dot(p,w)
r*w
o = uvtmp + r*w;

th = atan2(uvtmp(2)-o(2),uvtmp(1)-o(1)) + 0.02*pi()*linspace(-1,1,1000);
xc = o(1) + abs(r)*cos(th);
yc = o(2) + abs(r)*sin(th);
plot(xc, yc, 'k--');
axis image



r = abs(r);
d = norm(poly(1).xy(2,:) - o);
(d - r)/r
