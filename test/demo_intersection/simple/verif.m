clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');

xyz = importdata('curve_xyz.dat');
uv = cat(2, importdata('curve_uv1.dat'), importdata('curve_uv2.dat'));
cl = linspecerp(2);


q = round(0.33*size(xyz,1))


lbscale = 5e-1;
lboffset = 1e-3;
lbuv = -0.9*[1,1];

%%
figure;
hold on

for i = 1:2
    c = readCoeffs2(sprintf('c_%d.cheb',i));
    s = surf_chebyshev(c, 1);
    set(s, 'facecolor', cl(i,:));
    
    plot_local_basis( c, lbuv, lbscale, lboffset );
end

plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'k-', 'linewidth', 1);

plot3(xyz(q,1), xyz(q,2), xyz(q,3), 'k*');

axis image vis3d
view(47,25);
camproj('persp');
camlight(30,30);


%%
figure;

for i = 1:2
    subplot(1,2,i);
    hold on
    patch([-1,1,1,-1], [-1,-1,1,1], 'r', 'facecolor', cl(i,:));
    plot(uv(:,2*i-1), uv(:,2*i), 'k-', 'linewidth', 1);
    plot(uv(q,2*i-1), uv(q,2*i), 'k*');
    
    a = uv(q,2*i-1:2*i);
    v = uv(q+1,2*i-1:2*i) - a;
    v = 0.1*[-v(2), v(1)]/norm(v);
    plot(a(1), a(2), 'k*');
    quiver(a(1), a(2), v(1), v(2), 0, 'color', 'y', 'maxheadsize', 10);
    
    axis(repmat([-1,1],1,2));
    daspect([1,1,1])
end






return
%%
figure;
hold on

% for i = 1:2
%     c = readCoeffs2(sprintf('c_%d.cheb',i));
%     surf_chebyshev(c, 1);
% end

% points = importdata('intersection_points.dat');
% for ip = 1:size(points,1)
%     text(points(ip,1), points(ip,2), points(ip,3), num2str(ip), 'color', 'b', 'fontsize', 7);
% end

fid = fopen('intersection_curves.dat','r');
nc = str2num(fgetl(fid));
for ic = 1:nc
    fgetl(fid); % dummy
    fgetl(fid); % smooth
    fgetl(fid); % uvbox...
    fgetl(fid); % ...
    nsp = str2num(fgetl(fid)); % nsplit
    for isp = 1:nsp
        fgetl(fid); % isplit
    end
    for isp = 1:nsp-1
        fgetl(fid); % ...
    end
    np = str2num(fgetl(fid));
    uvxyz = zeros(np,7);
    for ip = 1:np
        uvxyz(ip,:) = str2num(fgetl(fid));
    end
%     plot3(uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), '-');
    a = uvxyz(1,5:7);
    v = uvxyz(np,5:7) - a;
    quiver3(a(1), a(2), a(3), v(1), v(2), v(3), 0, 'maxheadsize', 100);
    m = round(np/2);
    text(uvxyz(m,5), uvxyz(m,6), uvxyz(m,7), num2str(ic));
end
fclose(fid);

axis image vis3d
view(3);
camproj('persp');
camlight(30,30);
