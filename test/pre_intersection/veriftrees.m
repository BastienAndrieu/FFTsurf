clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');


cl = colorcet( 'I2', 'N', 4 );


for isurf = 1:2
    s(isurf).c = readCoeffs2(sprintf('c_%d.cheb',isurf));
    s(isurf).tree = importdata(sprintf('tree_%d.dat',isurf));
%     s(isurf).tree = importdata(sprintf('tree_%d_bezier.dat',isurf));
end


%% 2d
figure;
for isurf = 1:2
    subplot( 1, 2, isurf );
    hold on
    
    for j = 1:size(s(isurf).tree,1)
        rectangle( 'position', ...
            [ s(isurf).tree(j,[1,3]), s(isurf).tree(j,[2,4]) - s(isurf).tree(j,[1,3]) ] );
        
%         xm = 0.5*sum(s(isurf).tree(j,1:2));
%         ym = 0.5*sum(s(isurf).tree(j,3:4));
%         text( xm, ym, num2str(j), 'fontsize', 6, 'color', 'r' );
    end
    
%     axis image
    axis( repmat([-1,1],1,2) );
    daspect( [1,1,1] );
end



%% 3d
m = 100;

figure;
hold on

for isurf = 1:2
    
    si = surf_chebyshev( s(isurf).c, 1 );
    set( si, 'facecolor', cl(isurf,:) );
    
    if 1
        tree = s(isurf).tree;
        for i = 1:size(tree,1)
            uvi = [ linspace(tree(i,1),tree(i,2),m)' , tree(i,3)*ones(m,1) ;
                tree(i,2)*ones(m,1), linspace(tree(i,3),tree(i,4),m)' ;
                linspace(tree(i,2),tree(i,1),m)' , tree(i,4)*ones(m,1) ;
                tree(i,1)*ones(m,1), linspace(tree(i,4),tree(i,3),m)' ];
            
            c = ICT2unstr( s(isurf).c, uvi );
            plot3( c(:,1), c(:,2), c(:,3), 'k-' );
        end
    end
    
end

axis image vis3d
view(3)
camlight(30,30);
camproj('persp');








return

%% bug bernOBB2 (corrigé)


box(1).ctr = [0.188095659       9.03068285E-04   2.71757622E-03]';
box(1).rng = [-1.04787219     -0.763941348     -0.475386292]';
box(1).axe = [0.433332056     -0.176973507      0.883687556    
  0.884530842     -0.104389913     -0.454651445    
  0.172709331      0.978663921      0.111303031];


box(2).ctr = [0.235660687       1.49881635E-02   5.00041991E-02]';
box(2).rng = [-1.06843555     -0.653883040     -0.504000783]';
box(2).axe = [0.809905052     -0.547406077      0.210713878    
  0.263517648       1.86261367E-02 -0.964474738    
  0.524034560      0.836659789      0.159336567];

figure;
hold on
t = linspace(0,1,m)';

for isurf = 1:2
    h = plotOBB( box(isurf).ctr, box(isurf).rng, box(isurf).axe );
    set( h, 'facecolor', 'none', 'edgecolor', cl(isurf,:) );
    
    b = readCoeffs2(sprintf('bs_%d.bern',isurf));
    
    s = eval_Bezier_patch_tpgrid( b, t, t );
    surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
        'edgecolor','none','facecolor',cl(isurf,:),'diffusestrength',0.7, ...
        'facelighting','gouraud','specularstrength',0.15);
    
end

axis image vis3d
view(3)
camlight(30,30);
camproj('persp');
