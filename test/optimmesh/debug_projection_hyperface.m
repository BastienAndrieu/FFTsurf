clc; clear; close all

uv = [
    -0.99999998832876302      -0.34824230787771754
%     -0.95109041226696722       -4.9563153167452800E-002
%     -0.19982692784127659        1.0000000000000000
%     -0.90855108093640558      -0.87480251965961109
%     -0.82572867574825648      -0.95973446409775298
];
tng = [
    -13.136420291250859       -6.9352124739108234E-004
%     0.33042102111701488       -6.3406027484464120 
%     -1.0, 0.0
%     5.8718140794862750       -7.5289130120385606
%     6.1045728480497372       -2.6291021354910304
];
ab = [
    -2.2872891358716174E-009  -5.1100713655102104E-002
%     5.4145147668687166E-003  -6.3248970939102700E-002
%     -4.1238339000753799E-002   6.5815768343921000E-002
%     6.5708047934250002E-002  -2.2130114210398011E-002
];


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


norm(poly(1).xy(1,:) - poly(2).xy(end,:))



figure;
hold on

a = poly(1).xy(1,:);
v = poly(1).xy(2,:) - a;

quiver(a(1), a(2), v(1), v(2), ...
    0, 'b', 'maxheadsize', 1);
for i = 2
    plot(poly(i).xy(:,1), poly(i).xy(:,2), '.-');
    plot(poly(i).xy(1,1), poly(i).xy(1,2), 'o');
    for j = 1:ninter
        k = ipls(j,i);
        xy = (1.0 - lamb(j,i)) * poly(i).xy(k,:) + lamb(j,i) * poly(i).xy(k+1,:);
        plot(xy(1), xy(2), 'k*');
    end
end


l = 3*norm(v) / norm(tng);
quiver(uv(1), uv(2), l*tng(1), l*tng(2), 0, 'g', 'maxheadsize', 1);

axis image

tng(1)*v(2) - tng(2)*v(1)

[~,j] = max(lamb(:,1));
% j = 2;

omega = lamb(j,2);
fprintf('w = %f\n',omega);
k = ipls(j,2);
plot(poly(2).xy(k+[0,1],1), poly(2).xy(k+[0,1],2), 'go-');
% p = poly(2).xy(k+1,:) - poly(2).xy(k,:);
% if omega > 0.5
%     p = -p;
% end
% p = -p;
% tng = -tng;
p = ab;
quiver(uv(1), uv(2), p(1), p(2), 0, 'r', 'maxheadsize', 1);

w = [-tng(2), tng(1)] / norm(tng);

% r = 0.5 * sign(dot(tng,p)) * sum(p.^2) / dot(p,w)
r = 0.5 * dot(p,p) / dot(p,w);
r*w
o = uv + r*w;

th = atan2(uv(2)-o(2),uv(1)-o(1)) + 0.02*pi()*linspace(-1,1,1000);
xc = o(1) + abs(r)*cos(th);
yc = o(2) + abs(r)*sin(th);
plot(xc, yc, 'k-', 'markersize', 4);
axis(repmat(1.01*[-1,1],1,2));
daspect([1,1,1]);


r = abs(r);
d = norm(poly(1).xy(2,:) - o);
(d - r)/r

% ao = [
% -0.39556363413706924       0.86634801886978241
%     ];
% quiver(uv(1), uv(2), ao(1), ao(2), 0, 'r', 'maxheadsize', 1);
% plot(o(1), o(2), 'c*');

d = [
    -2.7766534497886131E-004   4.0980931548357316E-003
%     -2.3137220613229864E-003  -1.4437893682492575E-002
];
% quiver(uv(1), uv(2), d(1), d(2), 0, 'm', 'maxheadsize', 1);



w = 0.19414165594961694;
y = uv + w*ab;
% plot(y(1), y(2), 'm*');


p = [
    -0.91101555787263055      -0.87086925156361550
%     -0.82276049014832597      -0.96073413291137910
];
q = [
    -0.95136807761194608       -4.5465060012617069E-002  -6.0800762386375787E-002 -0.99962962259223687
%     -0.91131931163091606      -0.87120702856682297
%     -0.82295474840960725      -0.96091647233995370
];

% plot(p(1), p(2), 'r*');
% plot(q(1), q(2), 'k*');

return

pls = poly(2).xy(k+[0,1],:);
vec = pls(2,:) - pls(1,:);

vec = pls(2,:) - pls(1,:);
det([vec; poly(1).xy(2,:) - pls(1,:)])

w = 0.94023821046284350;%4.5172329619357687E-002;
wq = 0.94184700426482448;%1.8912562536821578E-003;


figure;
hold on

plot(poly(1).xy(2,1), poly(1).xy(2,2), '*');

xy = poly(2).xy;
xy = [xy(1:k,:); q; xy(k+1:end,:)];

plot(xy(:,1), xy(:,2), '.-');
% plot(pls(:,1), pls(:,2), 'go-');
% quiver(pls(1,1), pls(1,2), vec(1),vec(2), 0, 'maxheadsize', 1);

pls = [
    -0.94978941993715937      -0.80898675131569009     
 -0.91131931163091606      -0.87120702856682297
 ];

vec = pls(2,:) - pls(1,:);
det([vec; poly(1).xy(2,:) - pls(1,:)])



plot(pls(:,1), pls(:,2), 'ro-');

plot(pls(1,1), pls(1,2), pls(1,1) + w*vec(1), pls(1,2) + w*vec(2), '+');
plot(pls(1,1), pls(1,2), pls(1,1) + wq*vec(1), pls(1,2) + wq*vec(2), 's');

plot(q(1), q(2), 'o');

axis image
return

q = [
    -0.82284307326920514      -0.96096353230871656
];
figure;
hold on

plot(poly(1).xy(2,1), poly(1).xy(2,2), '*');

k = k + 1;
xy = [xy(1:k,:); q; xy(k+1:end,:)];

plot(xy(:,1), xy(:,2), '.-');
plot(pls(:,1), pls(:,2), 'go-');
quiver(pls(1,1), pls(1,2), vec(1),vec(2), 0, 'maxheadsize', 1);

axis image


return


figure;
hold on

plot(poly(2).xy(:,1), poly(2).xy(:,2), '.-');

fid = fopen('../polylines/xy3.dat','r');
np = str2num(fgetl(fid));
xy = zeros(np,2);
for j = 1:np
    xy(j,:) = str2num(fgetl(fid));
end
fclose(fid);

plot(xy(:,1), xy(:,2), '.-');

plot(poly(1).xy(2,1),poly(1).xy(2,2),'k*');

axis image

