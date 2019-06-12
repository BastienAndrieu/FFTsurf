clc; clear; close all

addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck8/bandrieu/Bureau/CYPRES/SphericalGeometry/');

%%
xyz = importdata('debuglonglat_xyz.dat');
A1 = importdata('debuglonglat_axes1.dat');
xyzp = importdata('debuglonglat_xyzproj.dat');
A2 = importdata('debuglonglat_axes2.dat');
xyp = importdata('debuglonglat_xyplane.dat');
obb = importdata('debuglonglat_obb.dat');

sqrt(sum(A1.^2,1)), sqrt(sum(A2.^2,1))


figure;
hold on

plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.');
plot3(xyzp(:,1), xyzp(:,2), xyzp(:,3), 'r*');

plot_triedre(A1, [0,0,0], 1);
plot_triedre(A2, [0,0,0], 1.5);

axis image vis3d
view(-160,20);
camproj('persp');


%%
figure;
hold on

plot(xyp(:,1), xyp(:,2), 'b.');

b = repmat(obb(1,:)',1,4) + obb(3:4,1:2)*(repmat(obb(2,:)',1,4).*[-1,1,1,-1;-1,-1,1,1]);
b = [b, b(:,1)];
plot(b(1,:), b(2,:), 'k--');

cl = [1,0,0;0,1,0];
for i = 1:2
    quiver(obb(1,1), obb(1,2), ...
        obb(2,i)*obb(3,i), obb(2,i)*obb(4,i), ...
        0, 'color', cl(i,:));
end

axis image


%%
uv = importdata('debuglonglat_uv.dat');

figure;
hold on

plot(uv(:,1), uv(:,2), 'b.');

axis image
% return

%%
uv = importdata('longlatpatch_uv.dat');
uv0 = importdata('/stck/bandrieu/Bureau/CYPRES/OffsetSurfaces/Huygens/corner_patch/bnd_uv.dat',' ',1);
uv0 = uv0.data;

c = readCoeffs2('longlatpatch_surf.cheb');

xyz0 = importdata('longlatpatch_xyz.dat',' ',1);
xyz0 = xyz0.data;
xyz = ICT2unstr(c, uv);


%%
figure;
hold on

surf_chebyshev(c, 1, 100, 5);
plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.');
plot3(xyz0(:,1), xyz0(:,2), xyz0(:,3), 'ro');

axis image vis3d
view(-160,20);
camlight(-30,30);
camproj('persp');

%%

% figure;
% hold on
% plot(uv(:,1), uv(:,2), 'b.');
% plot(uv0(:,1), uv0(:,2), 'ro');
% 
% axis(repmat(1.05*[-1,1], 1, 2));
% daspect([1,1,1]);



%%
N = size(c,1);
ctr = [0;1;1];
X = xyz0' - repmat(ctr, 1, size(xyz0,1));
[ C, uv0, Ts, Ps, solidangle ] = spherical_chebpatch_from_points( X, N );
C = 0.02 * C;
for i = 1:3
    C(1,1,i) = C(1,1,i) + ctr(i);
end

figure;
hold on

surf_chebyshev(C, 1, 100, 5);
plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.');
plot3(xyz0(:,1), xyz0(:,2), xyz0(:,3), 'ro');

axis image vis3d
view(-160,20);
camlight(-30,30);
camproj('persp');

uv0 = uv0';
figure;
hold on
plot(uv(:,1), uv(:,2), 'b.');
plot(uv0(:,1), uv0(:,2), 'ro');

axis(repmat(1.05*[-1,1], 1, 2));
daspect([1,1,1]);

