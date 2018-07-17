clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 16 );
[u,v] = tpgrid( cglpoints(-1,1,M), cglpoints(-1,1,N) );

%%
cls = linspecerp(2);
nisouv = 10;


%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            x = u;
            y = v;
            z = -0.6*(u.^2 + v.^2);
        case 2
            x = 0.5*u;
            y = 0.7*u.^2;
            z = v;
    end
    
    c = cat(3, fcht2d(x), fcht2d(y), fcht2d(z));
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,22) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


make_PN_coeff_matrices