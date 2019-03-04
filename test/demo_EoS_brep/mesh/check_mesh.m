clc; clear; close all

modevisu = 1;


if modevisu == 1
    ec = 'k';
    fl = 'flat';
else
    ec = 'none';
    fl = 'gouraud';
end

v = 6;
d = [-0.0000000000000000       -2.7289984793036839E-004  -0.0000000000000000];

% l = [2:10, 32:40, 162:180, 324:332, 432:440, 751:820];
l = 1;%2041,2072,10327,10358];%590,606];%602,637,2796,2831];
k = [10];%,121,86,120,54,127,41];

for i = 1%0:1
    xyz = importdata(sprintf('pos_%3.3d.dat',i));
    tri = importdata(sprintf('connect_%2.2d.dat',0));
    
    
    figure('units','normalized','position',[.1,.1,.6,.7])
    hold on
    
    trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
        'facecolor', 'w', 'edgecolor', ec, ...
        'facelighting', fl, 'specularstrength', 0);
    
    %     if i == 0
    %         plot3(xyz(v,1), xyz(v,2), xyz(v,3), 'bo');
    %         quiver3(xyz(v,1), xyz(v,2), xyz(v,3), d(1), d(2), d(3), 0, 'r');
    %     end
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'ro');
    
    axis image vis3d
    view(130,15);
    camlight(-30,30);
    camproj('persp');
    
    f = find(sum(tri == l(1),2))'
    for j = f
        b = (xyz(tri(j,1),:) + xyz(tri(j,2),:) + xyz(tri(j,3),:))/3;
        u = xyz(tri(j,2),:) - xyz(tri(j,1),:);
        v = xyz(tri(j,3),:) - xyz(tri(j,1),:);
        w = cross(u,v);
        quiver3(b(1), b(2), b(3), w(1), w(2), w(3), 'b');
    end
    
    trisurf(tri(50,:), xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'y');
    plot3(xyz(k,1), xyz(k,2), xyz(k,3), 'go');
    
    axis(xyz(l(1),[1,1,2,2,3,3]) + repmat(0.2*[-1,1],1,3));
    axis vis3d
    view(130,15);
    camlight(-30,30);
    camproj('persp');
end





R = 0.25;
for i = l
    ri = norm(xyz(i,:) - [xyz(i,1), 0, 0]);
    abs(ri - R)/R
end
