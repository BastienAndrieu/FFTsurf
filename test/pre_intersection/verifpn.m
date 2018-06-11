clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

cl = linspecerp( 2 );
[mc,nc] = deal(200,100);
u = linspace( 0, 2*pi(), mc )';



figure;
hold on
for isurf = 1:2
    pn = readCoeffs2( sprintf( 'pn_%d.cheb', isurf ) );
%     pn_eco = readCoeffs2( sprintf( 'pn_eco_%d.cheb', isurf ) );
    pn_eco = pn;
    
    %     [m,n,~] = size( pn );
    [p,q,~] = size( pn_eco);
    
    e = pn;
    e(1:p,1:q,:) = e(1:p,1:q,:) - pn_eco;
    
    %     surf_coeff_array( e );
    
    b = OBB_Chebyshev2( pn_eco );
    b.center
    b.extents
    b.axes
    
    
    [ep,fp] = plotOBBpyramid( b );
    
    set( ep, 'color', 'none' );
    set( fp, 'facecolor', cl(isurf,:), 'edgecolor', 'k' );
    
    
    xyz = obb2polyhedron( b );
    d = xyz ./ repmat( sqrt( sum(xyz.^2,2) ), 1, 3 );
    c = b.center' / norm(b.center);
    t = 0;
    for i = 1:size(d,1)
        ai = acos( min( 1, max( -1, dot(c,d(i,:)) ) ) );
        t = max( t, ai );
    end
    t/pi()
    
    [a,b] = completeBase( c' );
    R = [a,b,c'];
    
    r = tan(t)*linspace(0,1,nc);
    xc = cos(u)*r;
    yc = sin(u)*r;
    zc = repmat( linspace(0,1,nc), mc, 1 );
    
    xyzc = R * [ xc(:), yc(:), zc(:) ]';
    xyzc = reshape( xyzc', [mc, nc, 3] );
    
    surf( xyzc(:,:,1), xyzc(:,:,2), xyzc(:,:,3), ...
        'facecolor', min(1,1.2*cl(isurf,:)), 'edgecolor', 'none', ...
        'facelighting', 'gouraud', 'specularstrength', 0.1, 'ambientstrength', .6 );
end

view(3);
axis( repmat(1.05*[-1,1], 1, 3 ) );
daspect([1,1,1]);
axis vis3d
camlight(30,30);




