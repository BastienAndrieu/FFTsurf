clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

cl = ['b','r'];
cls = colorcet( 'I2', 'N', 2 );

mrk = ['+','x'];

clp = ['y','g'];



for isurf = 1:2
    %     S(isurf).bpn = readCoeffs2( sprintf('surfroot%d_pn.bern', isurf ) );
    S(isurf).bpn = readCoeffs2( sprintf('debugld_reg%d.bern', isurf ) );
    S(isurf).sep = importdata( sprintf( 'sep%d.dat', isurf ) );
end

[vec1, stat1] = separating_plane(S(1).sep,  S(2).sep)
[vec2, stat2] = separating_plane(S(1).sep, -S(2).sep)




rot = eye(3);
p = [
    1 0 0
    0 1 0
%     0.0000000000000000        4.3711389940630076E-008 -0.99999999999999900
%     -4.3711390001862412E-008   0.0000000000000000      -0.99999999999999900
    ]';
p = [vec1, vec2];

figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:2
    plot3( S(isurf).sep(:,1), S(isurf).sep(:,2), S(isurf).sep(:,3), ...
        '.', 'color', cl(isurf) );%, 'marker', mrk(isurf) );
end

for ip = 1:2
        [~,pla] = plot_separating_plane( p(:,ip) );
        set( pla, 'facecolor', clp(ip) );
end

cla = ['r','g','b'];
for i = 1:3
    quiver3( 0, 0, 0, rot(1,i), rot(2,i), rot(3,i), 0, 'color', cla(i) );
%     quiver3( 0, 0, 0, rot(i,1), rot(i,2), rot(i,3), 0, 'color', cla(i) );
end

for isurf = 1:2
    x = reshape( S(isurf).bpn, [], 3 );
    x = x ./ repmat( sqrt(sum(x.^2,2)), 1, 3 );
%     x = x * rot;
%     plot3( x(:,1), x(:,2), x(:,3), '.', 'color', cl(isurf), 'marker', mrk(isurf) );
end


axis( repmat( 1.05*[-1,1], 1, 3 ) );
daspect([1,1,1])
axis vis3d
view(3)
camlight(30,30);

%%
n = 200;
t = linspace(0,1,n)';

figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:2
    bs = S(isurf).bpn;
    if max(size(bs,1),size(bs,2)) < 2
        quiver3( 0, 0, 0, bs(1,1,1), bs(1,1,2), bs(1,1,3), 'color', cl(isurf) );
    else
        s = eval_Bezier_patch_tpgrid(bs, t, t);
        surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
            'edgecolor','none','facecolor',cls(isurf,:),'diffusestrength',0.7, ...
            'facelighting','gouraud','specularstrength',0.15);
    end
end

axis image vis3d
view(3)
camlight(30,30);
