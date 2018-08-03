clc; clear; close all

%%
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

cls = linspecerp(2);

nisouv = 5;


R = 1;
r = 0.2;
            
%%
M = 16;
[u,v] = tpgrid( cglpoints(1,-1,M) );
au = u;
av = v;
rot = 2*pi()*rand;

for isurf = 1:2
    switch isurf
        case 1 
            x = (R + r*cos(av)).*sin(au);
            y = r*sin(av);
            z = (R + r*cos(av)).*cos(au);
        case 2
            x = 2*(R+r)*u;
            y = 0*u;
            z = 2*(R+r)*v;
%             x = (R + r)*sin(av);
%             y = (R + r)*sin(au).*cos(av);
%             z = (R + r)*cos(au).*cos(av);
%             [x,y] = deal(cos(rot)*x + sin(rot)*y, -sin(rot)*x + cos(rot)*y);
    end
    S(isurf).C = cat( 3, fcht2d(x), fcht2d(y), fcht2d(z) );
    writeCoeffs2(S(isurf).C, sprintf('C%d.cheb', isurf));
end



%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    [si,isou,isov] = surf_chebyshev( S(isurf).C, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
%     if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
end


xlabel('x');
ylabel('y');
zlabel('z');

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);