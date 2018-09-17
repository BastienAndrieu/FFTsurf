clc; clear; close all

% mesh
tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
uv = importdata('uv.dat');
idstyp = importdata('idstyp.dat');

xyzs = importdata('xyz_smooth.dat');

% xyz = importdata('xyz_smooth.dat');
% idstyp = importdata('idstyp_smooth.dat');

% vertex displacement
dxyz = importdata('debug_dxyz.dat');

dxyz(491,:) = [
0.0000000000000000        1.7193342914554384E-003   8.4853348570061596E-004
];

figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

% quiver3(xyz(:,1), xyz(:,2), xyz(:,3), ...
%     dxyz(:,1), dxyz(:,2), dxyz(:,3), ...
%     0, 'color', 'r', 'maxheadsize', 1);

% l = 4243;
% l = 19408;
% l = 13887;%2045;
l = 491;
plot3( xyz(l,1), xyz(l,2), xyz(l,3), 'r*' );
quiver3(xyz(l,1), xyz(l,2), xyz(l,3), ...
    dxyz(l,1), dxyz(l,2), dxyz(l,3), ...
    0, 'color', 'r', 'maxheadsize', 1);
% k = 4261;
% k = [19409, 4261];
% k = [1936, 4261];
k = [22723,22724,22725,22726];
% plot3( xyz(k,1), xyz(k,2), xyz(k,3), 'g*' ); 


aabb = xyz(l,[1,1,2,2,3,3]) + 3e-2 * repmat([-1,1],1,3);



% axis image vis3d
axis(aabb); daspect([1,1,1]);
axis vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end


idstyp(l,:)
idstyp(k,:)


norm( xyz(l,:) - xyzs(l,:) )


% dxyz(l,:) * ( xyz(k,:) - repmat(xyz(l,:), length(k), 1) )'

uv(l,1:2)



figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyzs(:,1), xyzs(:,2), xyzs(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

plot3( xyzs(l,1), xyzs(l,2), xyzs(l,3), 'r*' );
% axis image vis3d
axis(aabb); daspect([1,1,1]);
axis vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end





return
% close all


x = [
    -8.0000000000000002E-002   5.8778525229246988E-002   8.0901699437494617E-002
    ];


d = [
    0.0000000000000000       -1.1799653577694766E-003   9.4490502347265124E-004
    ];

p = [
    -8.0000000000000002E-002   6.3263862061072512E-002   7.7444714197403855E-002
  -8.0000000000000002E-002   5.8778525229246988E-002   8.0901699437494617E-002
  ];



% J = [
%     -7.0685834705759912E-002  -1.7740149627076818E-014   0.0000000000000000     
%    0.0000000000000000        0.0000000000000000        7.5281646624104331E-002
%    0.16948314285097521        0.0000000000000000        1.7636199571454715E-002
%   -2.0816681711721685E-017   0.0000000000000000       -9.5281646624104321E-002
%    ];
% 
% 
% 
% J(1:2,:)' * [2.2981173103812625E-003 ; -1.4242238149147395E-002] - d'
% J(3:4,:)' * [-9.5846901115653651E-004 ; 1.1075326950213083E-002] - d'



figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

% for i = 1:2
%     quiver3(x(:,1), x(:,2), x(:,3), J(i,1), J(i,2), J(i,3), 'b');
% end

plot3(x(:,1), x(:,2), x(:,3), 'go');
quiver3(x(:,1), x(:,2), x(:,3), d(:,1), d(:,2), d(:,3), ...
    0, 'color', 'r', 'maxheadsize', 1);

plot3(p(:,1), p(:,2), p(:,3), 'y.-', 'markersize', 10);

axis(aabb); daspect([1,1,1]);
axis vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end







u = [
    -0.41993545977606123       0.85480037908316620
%     -0.27323954473512879        1.0000000000000000
    ]';

p = [
    -0.46213286263925812       0.83274407576811926     
 -0.41993545977606123       0.85480037908316620
%     -0.19366207318918069        1.0000000000000000     
%  -0.27323954473512879        1.0000000000000000
 ]';

du = [
    -8.1873061497877743E-003  -4.0406356461933496E-003
%     1.5257529859474579E-002   5.5436443983390216E-005
   ]';

tng = [
    -4.2740018954159398        7.9577471545958209     
  -2.0996772988800476        0.0000000000000000
  ];
tng = norm(du) * tng;


addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
% c = readCoeffs2('../../propergol/C_002.cheb');
% c = readCoeffs2('../../propergol/C_002.cheb');


sens = 1;

(-1)^sens * det([p(:,2) - p(:,1), u + du - p(:,1)])

figure;
hold on
if sens == 1
    quiver(p(1,2), p(2,2), p(1,1) - p(1,2), p(2,1) - p(2,2), 0, 'b', 'maxheadsize', 1)
else
    quiver(p(1,1), p(2,1), p(1,2) - p(1,1), p(2,2) - p(2,1), 0, 'b', 'maxheadsize', 1)
end
quiver(u(1), u(2), du(1), du(2), 0, 'r', 'maxheadsize', 1)

quiver(u(1), u(2), tng(1,sens), tng(2,sens), 0, 'g', 'maxheadsize', 1)

axis image

% t = linspace(0,2*pi(),1e3);
% r = 1 / 1.05;
% plot(r*cos(t), r*sin(t), 'k-');

return
% close all
% figure;
% hold on
% 
% c = readCoeffs2('../../propergol/C_002.cheb');
% surf_chebyshev(c, 1);
% 
% c = readCoeffs2('../../propergol/C_022.cheb');
% surf_chebyshev(c, 1);
% 
% plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'k*');
% 
% axis image vis3d
% view(3)
% camproj('persp');
% if 1
%     [azl,ell] = deal( 120, 60 );
%     [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
%     light( 'style', 'infinite', 'position', [xl,yl,zl] );
%     light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
% end

