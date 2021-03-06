clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');


plotfaces = 0;
plotdisplacement = 0;

instant = 61;

nsurf = 1;
cls = colorcet( 'I2', 'N', nsurf );
cls = cls(randperm(nsurf),:);
cls = CC( cls, 0.0, 0.8, 1.5 );

tri = importdata('../connect_01.dat');
xyz = importdata(sprintf('../pos_%3.3d.dat',instant));

if plotdisplacement
    d = importdata(sprintf('../pos_%3.3d.dat',instant+1));
    d = d - xyz;
end


figure('units','normalized', 'position', [.1 .1 .8 .8]);
hold on

if plotfaces
    for i = []%3,5,6,47,48]%1:nsurf
        if instant == 0
            c = readCoeffs2(sprintf('../../init/coef/c_%3.3d.cheb',i));
        else
            c = readCoeffs2(sprintf('instant_%3.3dc_%3.3d.cheb',instant,i));
        end
        si = surf_chebyshev(c, 1, 50);
        set(si, 'facecolor', cls(i,:) , 'specularstrength', 0);
    end
    trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'none', 'specularstrength', 0);
else
    if plotdisplacement
        magd = sqrt(sum(d.^2,2));
        trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), magd, ...
            'facecolor', 'interp', 'edgecolor', 'k');
    else
        trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'w', 'specularstrength', 0);
    end
end

if plotdisplacement
    quiver3(xyz(:,1), xyz(:,2), xyz(:,3), d(:,1), d(:,2), d(:,3), 0, 'r', 'maxheadsize', 1);
end

l = 23977;%[1943, 4379, 6705, 9031 ];%, 11202];%[1906,7293];%[16899];%
% plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'ro');


xyz = importdata(sprintf('../pos_%3.3d.dat',instant));

axis image vis3d
view(3)
% view(240,-50)
camproj('persp');
if plotfaces || ~plotdisplacement
    camlight(30,30);
end



addpath('/stck/bandrieu/Bureau/CYPRES/Regression3D/GPeyre/toolbox_graph/');
e = compute_edges(tri)';
le = sqrt(sum( (xyz(e(:,1),:) - xyz(e(:,2),:)).^2, 2));
figure;
hist(le)

mean(le)
min(le)
max(le)


% addpath('/stck/bandrieu/Bureau/CYPRES/Jouke2/cao/');
% [~,R] = sphereFit(xyz);
% R


addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/3d');
vol = volume_trisurf(tri, xyz);
vol
vol0 = 4.0*pi()*(0.15^3)/3.0;
abs(vol - vol0)/vol0
