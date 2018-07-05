clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

%%
uv = zeros(2,2,2);
xyz = zeros(3,2);
p = zeros(3,1);
fid = fopen( 'data.dat', 'r' );
fgetl(fid);%
for ipt = 1:2
    for isurf = 1:2
        uv(:,isurf,ipt) = str2num( fgetl(fid) )';
    end
end
fgetl(fid);%
for ipt = 1:2
    xyz(:,ipt) = str2num( fgetl(fid) )';
end
fgetl(fid);%
p = str2num( fgetl(fid) )';
fclose(fid);

ipt = 1;%randi(2);
uv = uv(:,:,ipt);
xyz = xyz(:,ipt);

for isurf = 1:2
    S(isurf).x = readCoeffs2( sprintf('c%d.cheb', isurf) );
    d = cheb_diff2( S(isurf).x );
    S(isurf).xu = d(:,:,:,1);
    S(isurf).xv = d(:,:,:,2);    
end
%%
uvc = [1,1];
scln = 0.3;

cls = colorcet( 'I2', 'N', 2 );


figure;
hold on
for isurf = 1:2   
    si = surf_chebyshev( S(isurf).x, 1 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
    
    x = ICT2unstr( S(isurf).x, uvc );
    xu = ICT2unstr( S(isurf).xu, uvc );
    xv = ICT2unstr( S(isurf).xv, uvc );
    n = cross( xu, xv );
    n = scln * n / norm( n );
    
    quiver3( x(1), x(2), x(3), n(1), n(2), n(3), 0, 'color', 0.5*cls(isurf,:) );
end

plot3( xyz(1), xyz(2), xyz(3), 'k*' );
quiver3( xyz(1), xyz(2), xyz(3), p(1), p(2), p(3), 'r' );


axis image vis3d
daspect( [1,1,1] );
view(3)
camproj('persp');
camlight(30,30);




