clc; clear; close all

addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck8/bandrieu/Bureau/CYPRES/Intersections');

cl = colorcet( 'I2', 'N', 2 );
cl = cl(randperm(size(cl,1)),:);
cl = CC( cl, 0.0, 0.8, 1.5 );


uv = [
    0.99999999999999689       0.98009925619580018
    -0.99999999999999689       0.98009925619580018
    ];



xyz = zeros(1,3);
for i = 1:2
    c = readCoeffs2(sprintf('debug_diffgeominter_surf%d.cheb',i));
    d = cheb_diff2(c);
    S(i).D = d;
    d2u = cheb_diff2(d(:,:,:,1));
    d2v = cheb_diff2(d(:,:,:,2));
    S(i).D2 = cat(4, d2u, d2v(:,:,:,1));
    
    si = surf_chebyshev(c,1,200);
    set(si, 'facecolor', cl(i,:), 'specularstrength', 0);
    
    xyzi = ICT2unstr(c, uv(i,:));
    xyz = xyz + 0.5 * xyzi;
end

plot3(xyz(1), xyz(2), xyz(3), 'k*');

axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );

[ stat, duv_ds, dxyz_ds ] = diffgeom_intersection_curve( S, uv )

for i = 1:size(dxyz_ds,2)
    quiver3(xyz(1), xyz(2), xyz(3), dxyz_ds(1,i), dxyz_ds(2,i), dxyz_ds(3,i), ...
        'k' );
end

