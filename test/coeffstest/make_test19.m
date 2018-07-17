clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');

%%
cls = linspecerp(2);
nisouv = 10;

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            f = 5;
            x = [
                -1,-1,-1,-1;
                1 , 1, 1, 1;
                -1,-1,-1,-1;
                1 , 1, 1, 1];
            x(2:3,:) = f*x(2:3,:);            
            y = x';
            
            v = linspace(-1,1,4);
            z = 0.05 * ( repmat(v',1,4) + repmat(v,4,1) );
            
        case 2
            x = repmat( scale2*linspace(-1,1,4)', 1, 4 );
            y = x';
            fz = 0.6;
            z = 0.6 * [
                1    -5     5    -1
                -5    25   -25     5
                5   -25    25    -5
                -1     5    -5     1];

    end
    
    b = cat( 3, x, y, z );
    c = 0*b;
    A = C2Bmatrix(4);
    At = A';
    for i = 1:3
        c(:,:,i) = A \ b(:,:,i) / At;
    end
    
%     C = cat( 3, xc, yc, zc );
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,19) );
    
    %     B = chebyshev2bezier_2( C );
    %     surf( B(:,:,1), B(:,:,2), B(:,:,3), 'facecolor', 'none', 'edgecolor', 0.2*cls(isurf,:) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


make_PN_coeff_matrices