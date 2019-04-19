clc; clear; close all


fid = fopen('../debug/faces_new.dat');
nf = 0;
while ~feof(fid)
    fgetl(fid);
    nin = str2num(fgetl(fid));
    for iin = 1:nin
        fgetl(fid);
    end
    nf = nf + 1;
end

cl = colorcet( 'I2', 'N', nf );
cl = CC( cl, 0.0, 0.8, 1.5 );
% cl = cl(randperm(size(cl,1)),:);
% if 0
%     save('/d/bandrieu/GitHub/These/memoire/figures/code/demo_EoS_brep_palette.dat',...
%         'cl','-ascii');
% end
% cl = importdata('/d/bandrieu/GitHub/These/memoire/figures/code/demo_EoS_brep_palette.dat');
cliso = 0.1*cl;

modevisu = 2;


if modevisu == 1
    ec = 'k';
    fl = 'flat';
else
    ec = 'none';
    fl = 'gouraud';
end

%%
figure;
hold on

for iface = 1:nf
    c = readCoeffs2(sprintf('c_%3.3d.cheb',iface));
            [si, lu, lv] = surf_chebyshev(c, 1, 100, 9);
            set(si, 'facecolor', cl(iface,:));
            set([lu;lv], 'color', cliso(iface,:));
end

axis image vis3d
view(130,15);
camlight(-30,30);
camproj('persp');

return

%%
figure;
hold on

for iface = 1:nf
    xyz = importdata(sprintf('xyz_%3.3d.dat',iface));
    tri = importdata(sprintf('tri_%3.3d.dat',iface));
    trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
        'facecolor', cl(iface,:), 'edgecolor', ec, ...
        'facelighting', fl, 'specularstrength', 0);
end

axis image vis3d
view(130,15);
camlight(-30,30);
camproj('persp');

%%
ratio = 1.5;
row = sqrt(nf/ratio);
row = max(1, round(row));
col = ceil(nf/row);

figure;
for iface = 1:nf
    uv = importdata(sprintf('uv_%3.3d.dat',iface));
    tri = importdata(sprintf('tri_%3.3d.dat',iface));
    
    subplot(row,col,iface);
    trisurf(tri, uv(:,1), uv(:,2), 0*uv(:,1), 'facecolor', cl(iface,:));
    axis(repmat(1.05*[-1,1], 1, 2));
    daspect([1,1,1]);
end

