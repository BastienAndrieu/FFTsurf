clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

surftag = importdata('surftag.dat');
nsurf = surftag(1);
cls = colorcet( 'I2', 'N', nsurf );
cls = CC( cls, 0.0, 0.8, 1.5 );

for isurf = 1:nsurf
    face(isurf).c = readCoeffs2(sprintf('coef/c_%3.3d.cheb',isurf));
end

fid = fopen('tangent_curves.dat','r');
nc = str2num(fgetl(fid));
for ic = 1:nc
    curves(ic).surf = str2num(fgetl(fid));
    np = str2num(fgetl(fid));
    curves(ic).xyz = zeros(np,3);
    for i = 1:np
        curves(ic).xyz(i,:) = str2num(fgetl(fid));
    end
    curves(ic).uv = zeros(np,4);
    for i = 1:np
        curves(ic).uv(i,:) = str2num(fgetl(fid));
    end
    curves(ic).uv = reshape(curves(ic).uv, [np,2,2]);
end
fclose(fid);


figure;
hold on

for isurf = 1:nsurf
    si = surf_chebyshev(face(isurf).c, 1);
    set(si, 'facecolor', cls(isurf,:) , 'specularstrength', 0);
end

for ic = 1:nc
    plot3(curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
        'k.-');
end

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
