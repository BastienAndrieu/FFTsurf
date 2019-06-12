clc; clear; close all

tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
v2h = importdata('mv2h.dat');
twin = importdata('mtwin.dat');

figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

% l = find(min(twin,[],2) < 1);
% l = find(v2h < 1);
% plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'r.');

l = [];
for i = 1:size(twin,1)
    b = sum(xyz(tri(i,:),:),1) / 3.0;
    for j = 1:3
        if twin(i,2*j) < 1
            l = [l; [j,i]];
            e = tri(i,[j, 1+mod(j,3)])
            a = xyz(e(1),:);
            v = xyz(e(2),:) - a;
            m = a + 0.5*v;
            a = a + 0.05 * (b - m) + 0.05 * v;
            v = 0.9*v;
            quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, 'r', 'maxheadsize', 1);
        end
    end
end
l


axis image vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end

% return

nt = size(tri,1);
nv = size(xyz,1);

j = 40;%randi(nv);
% j = l(randi(length(l)));

h = v2h(j,:);
kv = [];
kf = [];
while 1
    % face
    kf = [kf, h(2)];
    % destination
    dest = tri(h(2), 1+mod(h(1),3));
    kv = [kv, dest];
    % previous (ingoing)
    h(1) = 1 + mod(h(1)+3-2,3);
    % twin (outgoing)
    h = twin(h(2),2*h(1)-[1,0]);
    if all(h == v2h(j,:)) || any(h == 0)
        break
    end
end

l = [j,kv];
aabb = zeros(1,6);
aabb(1:2:5) = min(xyz(l,:),[],1);
aabb(2:2:6) = max(xyz(l,:),[],1);
rng = aabb(2:2:6) - aabb(1:2:5);
rng = max(rng, 0.1*max(rng));
aabb(1:2:5) = aabb(1:2:5) - rng;
aabb(2:2:6) = aabb(2:2:6) + rng;


nor = zeros(1,3);
for f = kf
    nor = cross(xyz(tri(f,2),:) - xyz(tri(f,1),:), xyz(tri(f,3),:) - xyz(tri(f,1),:));
end
nor = 0.4 * norm(rng) * nor / norm(nor);



figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

plot3(xyz(j,1), xyz(j,2), xyz(j,3), 'r.');
for i = 1:length(kv)
    a = xyz(kv(i),:);
    b = xyz(kv(1+mod(i,length(kv))),:);
    v = b - a;
    quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, 'r', 'maxheadsize', 1);
end

quiver3(xyz(j,1), xyz(j,2), xyz(j,3), nor(1), nor(2), nor(3), 0, 'b', 'maxheadsize', 1);

trisurf(tri(kf,:), xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'y', 'specularstrength', 0, 'edgecolor', 'none');
for i = 1:length(kf)
    fi = kf(i);
    
    a = sum(xyz(tri(kf(i),:),:),1) / 3.0;
    b = sum(xyz(tri(kf(1+mod(i,length(kf))),:),:),1) / 3.0;
    v = b - a;
    quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, 'g', 'maxheadsize', 1);
end

% axis image vis3d
axis(aabb);
daspect([1,1,1]);
axis vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end