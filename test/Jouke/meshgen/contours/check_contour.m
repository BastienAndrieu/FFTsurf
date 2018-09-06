clc; clear; close all

iface = 55;



fid = fopen(sprintf('uv_%3.3d.dat',iface), 'r');
nd = str2num(fgetl(fid));
n = nd(1);
d = nd(2);
x = zeros(n,d);
for i = 1:n
    x(i,:) = str2num(fgetl(fid));
end
fclose(fid);


fid = fopen(sprintf('edges_%3.3d.dat',iface), 'r');
ne = str2num(fgetl(fid));
e = zeros(ne,2);
for i = 1:ne
    e(i,:) = str2num(fgetl(fid));
end
fclose(fid);





figure;
hold on

plot(x(:,1), x(:,2), '.');
% for i = 1:n
%     text(x(i,1), x(i,2), num2str(i), 'color', 'r');
% end


% 
axis image
% 
for i = 1:ne
    plot(x(e(i,:),1), x(e(i,:),2), 'k-');
end


% return


uv = importdata(sprintf('../uv_%3.3d.dat',iface));
tri = importdata(sprintf('../tri_%3.3d.dat',iface));


figure;
hold on
triplot(tri, uv(:,1), uv(:,2) );
for i = 1:ne
    plot(x(e(i,:),1), x(e(i,:),2), 'r.-');
end
axis image



%%
c = readCoeffs2(sprintf('../coeffs/c_%3.3d.cheb',iface));

figure;
hold on

surf_chebyshev(c, 1, 200, 11);

xyz = ICT2unstr(c, uv);
trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'none');

axis image vis3d
view(3)
camlight(30,30);
