clc; clear; %close all

verts = importdata('verts.dat');
nv = size(verts,1);
edges = importdata('edges.dat'); % face, orig, prev(2), next(2)

%%
fid = fopen( 'intersection_curves.dat', 'r' );
ne = str2num( fgetl( fid ) );
for ic = 1:ne
    curves(ic).dummy = str2num( fgetl( fid ) );
    curves(ic).smooth = str2num( fgetl( fid ) );
    curves(ic).uvbox = zeros( 2, 4 );
    for isurf = 1:2
        curves(ic).uvbox(isurf,:) = str2num( fgetl( fid ) );
    end
    curves(ic).nsplit = str2num( fgetl( fid ) );
    curves(ic).isplit = zeros( curves(ic).nsplit, 2 );
    for ip = 1:curves(ic).nsplit
        curves(ic).isplit(ip,:) = str2num( fgetl( fid ) );
    end
    curves(ic).class = zeros( curves(ic).nsplit-1, 1 );
    for is = 1:curves(ic).nsplit-1
        curves(ic).class(is) = str2num( fgetl( fid ) );
    end
    n = str2num( fgetl( fid ) );
    curves(ic).np = n;
    curves(ic).uv = zeros( n, 2, 2 );
    curves(ic).xyz = zeros( n, 3 );
    for ip = 1:n
        line = str2num( fgetl( fid ) );
        curves(ic).uv(ip,:,1) = line(1:2);
        curves(ic).uv(ip,:,2) = line(3:4);
        curves(ic).xyz(ip,:) = line(5:7);
    end
end
fclose( fid );

%%
fid = fopen('edges_xyz.dat','r');
ne = str2num(fgetl(fid));
for ie = 1:ne
    np = str2num(fgetl(fid));
    edgexyz(ie).xyz = zeros(np,3);
    for ip = 1:np
        edgexyz(ie).xyz(ip,:) = str2num(fgetl(fid));
    end
end
fclose(fid);

%%
fid = fopen('faces.dat');
nf = 0;
while ~feof(fid)
    fgetl(fid);
    nin = str2num(fgetl(fid));
    for iin = 1:nin
        fgetl(fid);
    end
    nf = nf + 1;
end

%%
cl = colorcet( 'I2', 'N', nf );
cl = CC( cl, 0.0, 0.8, 1.5 );

modevisu = 2;
if modevisu == 1
    ec = 'k';
    fl = 'flat';
else
    ec = 'none';
    fl = 'gouraud';
end


figure('color', 'w');
hold on

for iface = 1:nf
    xyz = importdata(sprintf('../brepmesh/xyz_%3.3d.dat',iface));
    tri = importdata(sprintf('../brepmesh/tri_%3.3d.dat',iface));
    trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
        'facecolor', cl(iface,:), 'edgecolor', ec, ...
        'facelighting', fl, 'specularstrength', 0, ...
        'facealpha', 0.5, 'ambientstrength', 0.3);
end

for iedge = 1:ne
    %     xyz = curves(iedge).xyz;
    xyz = edgexyz(iedge).xyz;
    plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'k-');
    i = round(0.5*size(xyz,1));
    a = xyz(i,:);
    v = xyz(i+1,:) - a;
    v = 0.1*v/norm(v);
    quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, 'k', 'maxheadsize', 10);
    text(a(1), a(2), a(3), num2str(iedge));
end

for ivert = 1:nv
    xyz = verts(ivert,1:3);
    text(xyz(1), xyz(2), xyz(3), num2str(ivert), 'color', 'r', 'fontweight', 'bold');
end

axis image vis3d off
view(130,15);
camlight(-30,30);
camproj('persp');
