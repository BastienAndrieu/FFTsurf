clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

for isurf = 1:2
    ichild = randi(4);
    b = readCoeffs2(sprintf('child_%d_%d_pn.bern',isurf,ichild));
    
    gm(isurf).b = reshape( b ./ repmat(sqrt(sum(b.^2,3)),[1,1,3]), [], 3 );
    
end

cl = ['r','b'];

figure;
hold on

for isurf = 1:2
    plot3( gm(isurf).b(:,1), gm(isurf).b(:,2), gm(isurf).b(:,3), '.', 'color', cl(isurf) );
end

axis( repmat( 1.05*[-1,1], 1, 2 ) );
daspect([1,1,1]);
view(3)