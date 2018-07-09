clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');


xc = [ 0, 0, 0, 0;
    1, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0];
yc = [ 0, 1, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0];
scale2 = 0.9;

%%
cls = linspecerp(2);
nisouv = 10;

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            zc = [ 0, 0, 0, 0;
                0, 0, 0, 0;
                0, 0, 0, 0;
                0, 0, 0, 0];
        case 2
            zc = [ 0.01, 0, 0, 0;
                0, 0, 0, 0;
                0, 0, 0, 0;
                0, 0, 0, -0.2];
            xc = scale2*xc;
            yc = scale2*yc;
    end
    
    C = cat( 3, xc, yc, zc );
    
    [si,isou,isov] = surf_chebyshev( C, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( C, sprintf('C%d_test%2.2d.txt',isurf,18) );
    
%     B = chebyshev2bezier_2( C );
%     surf( B(:,:,1), B(:,:,2), B(:,:,3), 'facecolor', 'none', 'edgecolor', 0.2*cls(isurf,:) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


make_PN_coeff_matrices