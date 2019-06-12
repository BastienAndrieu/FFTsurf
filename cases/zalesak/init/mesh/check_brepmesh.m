clc; clear; close all

tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
idstyp = importdata('idstyp.dat');


nt = size(tri,1);
nv = size(xyz,1);

nverts = 0;
nedges = 0;
nfaces = 0;
for i = 1:nv
    switch idstyp(i,2)
        case 2
            nfaces = max(nfaces, idstyp(i,1));
        case 1
            nedges = max(nedges, idstyp(i,1));
        case 0
            nverts = max(nverts, idstyp(i,1));
    end
end

fid = fopen('paths.dat','r');
npths = str2num(fgetl(fid));
for i = 1:npths
    pth(i).hyperedge = str2num(fgetl(fid));
    pth(i).nv = str2num(fgetl(fid));
    pth(i).verts = zeros(1,pth(i).nv);
    pth(i).s = zeros(pth(i).nv,1);
    for j = 1:pth(i).nv
        l = str2num(fgetl(fid));
        pth(i).verts(j) = round(l(1));
        pth(i).s(j) = l(2);
    end
end
fclose(fid);





% nfaces = 179;%
cl = colorcet( 'I2', 'N', nfaces );
cl = cl(randperm(size(cl,1)),:);
% cl = CC( cl, 0.0, 0.8, 1.5 );

clp = 0.8 * colorcet( 'R3', 'N', max(2,npths) );

figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');


% trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ener, ...
%     'specularstrength', 0, 'edgecolor', 'none');

l = 4060:4093;%1936:2044;%1:65;
plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.', 'color', 'g', 'markersize', 6);

if 0
    %     for iface = 1:nfaces
    %         l = find(idstyp(:,1) == iface & idstyp(:,2) == 2);
    %         plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.', 'color', cl(iface,:), 'markersize', 6);
    %     end
    
    l = find(idstyp(:,2) == 1);
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.', 'color', 'g', 'markersize', 10);
    
    l = find(idstyp(:,2) == 0);
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'o', 'color', 'k', 'markersize', 6);
end

if 0
    for i = 1:npths
        l = pth(i).verts;
        %     plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'o', 'color', clp(i,:), 'linewidth', 2);
        plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'o-', 'color', clp(i,:), 'linewidth', 2);
        %     plot3(xyz(l(1),1), xyz(l(1),2), xyz(l(1),3), '*', 'color', 'r', 'markersize', 10);
        %         plot3(xyz(l(2),1), xyz(l(2),2), xyz(l(2),3), '*', 'color', 'y', 'markersize', 10);
        %     plot3(xyz(l(end),1), xyz(l(end),2), xyz(l(end),3), 'o', 'color', clp(i,:), 'markersize', 8);
    end
end

axis image vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end


%%
figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

% l = randi(nt);
l = find(sum(tri == 4060,2));
l = l(randi(length(l)));
f = tri(l,:);
patch(xyz(f,1), xyz(f,2), xyz(f,3), 'y');

for k = 1:3
    i = f(k);
    j = f(1+mod(k,3));
    a = xyz(i,:);
    v = xyz(j,:) - a;
    quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, 'r', 'maxheadsize', 1);
end

axis(a([1,1,2,2,3,3]) + 0.1*repmat([-1,1],1,3));
daspect([1,1,1])
axis vis3d
view(3)
camproj('persp');
camlight(30,30);