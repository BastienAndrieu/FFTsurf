clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

uv = [
    1.0000000000000000      -0.21830194096289457       0.17106892505234828       0.20576203110265048
    -0.68099312119192201       0.15832606383026293       -1.0000000000000000      -0.89074736846327596];

xyz = [
    -1.3768309561549987       -3.4693106050302201E-002   7.9324719935400884E-002
    0.89074736846330926      -0.10925263153662537       0.12154675604187920];

p = xyz(2,:) - xyz(1,:);

fid = fopen( 'data.dat', 'w' );
fprintf( fid, 'uv_endpoints\n');
for ipt = 1:2
    for isurf = 1:2
        fprintf( fid, '%16.15e %16.15e\n', uv(ipt,2*isurf-1), uv(ipt,2*isurf) );
    end
end
fprintf( fid, 'xyz_endpoints\n');
for ipt = 1:2
    fprintf( fid, '%16.15e %16.15e %16.15e\n', xyz(ipt,1), xyz(ipt,2), xyz(ipt,3) );
end
fprintf( fid, 'param_vector\n');
    fprintf( fid, '%16.15e %16.15e %16.15e\n', p(1), p(2), p(3) );
fclose(fid);


cls = colorcet( 'I2', 'N', 2 );

figure;
hold on
for isurf = 1:2
    c = readCoeffs2( sprintf( ...
        '/stck/bandrieu/Bureau/coeffstest/C%d_test01.txt', isurf) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
    
    writeCoeffs2( c, sprintf('c%d.cheb', isurf) );
end

plot3( xyz(:,1), xyz(:,2), xyz(:,3), 'k*' );

for i = 1:2
    quiver3( xyz(i,1), xyz(i,2), xyz(i,3), p(1), p(2), p(3), 'r' );
end

axis image vis3d
daspect( [1,1,1] );
view(3)
camproj('persp');
camlight(30,30);