clc; clear; close all

tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
idstyp = importdata('idstyp.dat');

ener = importdata('ener.dat');
grad = importdata('grad.dat')';
hess = importdata('hess.dat');

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

% l = 1936:2044;%1:65;
% plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.', 'color', 'g', 'markersize', 6);

if 1
%     for iface = 1:nfaces
%         l = find(idstyp(:,1) == iface & idstyp(:,2) == 2);
%         plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.', 'color', cl(iface,:), 'markersize', 6);
%     end
    
    l = find(idstyp(:,2) == 1);
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.', 'color', 'g', 'markersize', 10);
    
    l = find(idstyp(:,2) == 0);
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'o', 'color', 'k', 'markersize', 6);
end

for i = 1%:npths
    l = pth(i).verts;
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'o', 'color', clp(i,:), 'linewidth', 2);
%     plot3(xyz(l,1), xyz(l,2), xyz(l,3), '-', 'color', clp(i,:), 'linewidth', 2);
    plot3(xyz(l(1),1), xyz(l(1),2), xyz(l(1),3), '*', 'color', 'r', 'markersize', 10);
        plot3(xyz(l(2),1), xyz(l(2),2), xyz(l(2),3), '*', 'color', 'y', 'markersize', 10);

%     plot3(xyz(l(end),1), xyz(l(end),2), xyz(l(end),3), 'o', 'color', clp(i,:), 'markersize', 8);
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

return

%%

figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

for i = 1:npths
    l = pth(i).verts;
    plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.-', 'color', clp(i,:), 'linewidth', 1);
    
    m = 1;%round(length(l)/2);
    a = xyz(l(m),:);
    b = xyz(l(m+1),:);
    v = b - a;
    v = 0.1 * v / norm(v);
    
    quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, ...
        'color', clp(i,:), 'linewidth', 1.5, 'MaxHeadSize', 100 );
end
    
axis image vis3d
view(3)
camproj('persp');

return


%%
for i = 1:nv
    if idstyp(i,2) == 2
        d = - reshape(hess(i,:),3,3) \ grad(:,i);
        xyz(i,:) = xyz(i,:) + d';
    end
end
figure('units','normalized','position',[.15 .15 .7 .7]);
hold on
trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');
axis image vis3d
view(3)
camproj('persp');
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );

return

fid = fopen('../../brep/hyperfaces.dat','r');
nhf = str2num(fgetl(fid));
for j = 1:nhf
    hyperface(j).nf = str2num(fgetl(fid));
    hyperface(j).faces = str2num(fgetl(fid));
end
fclose(fid);

fid = fopen('../../brep/hyperedges.dat', 'r');
nhe = str2num(fgetl(fid));
for j = 1:nhe
    hyperedge(j).ne = str2num(fgetl(fid));
    hyperedge(j).verts = str2num(fgetl(fid));
    hyperedge(j).hyperfaces = str2num(fgetl(fid));
    for k = 1:hyperedge(j).ne
        hyperedge(j).edges(k,:) = str2num(fgetl(fid));
    end
end
fclose(fid);


f2hf = zeros(1,nfaces);
for j = 1:nhf
    for k = hyperface(j).faces
        f2hf(k) = j;
    end
end

e2he = zeros(1,nedges);
for j = 1:nhe
    for k = hyperedge(j).edges
        e2he(k) = j;
    end
end

refv = zeros(nv,1);
for i = 1:nv
    switch idstyp(i,2)
        case 2
            refv(i) = f2hf(idstyp(i,1));
        case 1
            refv(i) = e2he(idstyp(i,1));
        case 0
            refv(i) = idstyp(i,1);
    end
end

% refe = zeros(ne,1);




return
addpath('/stck8/bandrieu/Bureau/CYPRES/Regression3D/GPeyre/toolbox_graph/');
edg = compute_edges(tri');
ne = size(edg,2);

fprintf('Euler-Poincaré = %d\n', nv - ne + nt);

addpath('/stck8/bandrieu/Bureau/INRIA/Matlab/');
inria_mesh( 'brepmesh.mesh', ...
    xyz', tri', edg, refv, zeros(nt,1), zeros(ne,1) );