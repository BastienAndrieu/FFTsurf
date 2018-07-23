clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');

%%
cls = linspecerp(2);
nisouv = 10;

pair = [15,2];

scale = 15;

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    
    c = readCoeffs2(sprintf(...
        '/stck/bandrieu/Bureau/FFTsurf/cases/Jouke2/reference/Patchs/C_%2.2d.cheb',pair(isurf)));
    c = c * scale;
    if isurf == 1
        t = zeros(3,1);%c(1,1,:) / scale;
    end
    for idim = 1:3
        c(:,:,idim) = c(:,:,idim) - t(idim);
    end
    
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,25) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');