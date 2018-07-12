clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

cl = colorcet( 'I2', 'N', 2 );

uvbox = [
    -1.0000000000000000       0.21100575942548269       -1.0000000000000000      -0.18439541145683380
    -0.43041103318309587       0.13917793363380829      -0.58740087525063678      -0.17480175050127350
    ];


[icurv, ivar, ival] = deal(2,2,1);
jvar = 1 + mod(ivar,2);

% tuvxyz = [
%     -0.65467850697217500      -0.22256436409996938      -0.64470027716729394      -0.81895398933095631      -0.39918135451747871        2.3719879497097356E-002; % 2,2,1
%     -4.7996571522888734E-002 -0.31807366295874689      -0.43758219926467290      -0.50806480844753121      -0.55982811631006035        5.0220268170469914E-002; % 2,1,1
%     ];

uvxyz = [
    -0.37923765875335397      -0.21331415834057010      -0.48514793768139386      -0.17480175050127350      -0.17506027735847401      -0.64591440397182875        3.2407933135260830E-002
    ];

t = uvxyz(2*(icurv-1)+jvar);
t = -1 + 2 * (t - uvbox(icurv,2*jvar-1)) / (uvbox(icurv,2*jvar) - uvbox(icurv,2*jvar-1))
erruv = uvxyz(2*(icurv-1)+ivar) - uvbox(icurv,2*(ivar-1)+ival)


figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    c = chgvar2( c, reshape(uvbox(isurf,:),2,2) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    
    if isurf == icurv
        c = border_parameterization( c, ivar, ival );
        g = chebval1( c, linspace(-1,1,100)' );
        plot3( g(:,1), g(:,2), g(:,3), 'k-' );
    end
end




plot3( uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), 'k.', 'markersize', 8 );

plot3( -0.50806480844753121,-0.55982811631006035,5.0220268170469914E-002,'b*');
plot3( -0.81895398933095631,-0.39918135451747871,2.3719879497097356E-002,'r*');

axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );

