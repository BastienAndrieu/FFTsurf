clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

%%
cls = linspecerp(2);
nisouv = 10;

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    
    c = readCoeffs2(sprintf('C%d_test01.txt',isurf));
    x = ifcht2d(c(:,:,1));
    y = ifcht2d(c(:,:,2));
    z = ifcht2d(c(:,:,3));
    y = -y;
    c = cat(3, fcht2d(x), fcht2d(y), fcht2d(z));
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,21) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


make_PN_coeff_matrices