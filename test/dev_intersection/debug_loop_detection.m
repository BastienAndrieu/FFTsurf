clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');

cl = ['b','r'];
cls = colorcet( 'I2', 'N', 2 );
clp = ['y','g'];

for isurf = 1:2
    S(isurf).bpn = readCoeffs2( sprintf('surfroot%d_pn.bern', isurf ) );
    S(isurf).sep = importdata( sprintf( 'sep%d.dat', isurf ) );
end


p = [0.0000000000000000      -0.81902447538258794       0.57375858052344342;
    0.0000000000000000       0.86135953765087536       0.50799581385860859]';


figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:2
    plot3( S(isurf).sep(:,1), S(isurf).sep(:,2), S(isurf).sep(:,3), ...
        '.', 'color', cl(isurf) );
end

for ip = 1:2
    [~,pla] = plot_separating_plane( p(:,ip) );
    set( pla, 'facecolor', clp(ip) );
end

axis( repmat( 1.05*[-1,1], 1, 3 ) );
axis vis3d
view(3)
camlight(30,30);