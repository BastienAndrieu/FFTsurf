clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
addpath('/stck/bandrieu/Bureau/CYPRES/SphericalGeometry/');


cls = linspecerp(2);


% p1 = [0.55279632227442066       0.0000000000000000        0.83331640214259239]';
% p2 = [0.0000000000000000        0.0000000000000000        1.0000000000000000]';
% uvb = [0.62464114110311919        1.0000000000000000      -0.70211437915684782       0.14894281042157609     
%   0.33558402550468880        1.0000000000000000       -1.0000000000000000      -0.94656563208920530];

% p1 = [0.91013776496971577        0.0000000000000000      -0.41430574311241503]';
% p2 = [0.0000000000000000        0.0000000000000000        1.0000000000000000]';
% uvb = [-0.37535885889688086       0.24928228220623827       0.14894281042157609        1.0000000000000000
%     0.33558402550468880        1.0000000000000000        5.3434367910794645E-002   1.0000000000000000];
    

p1 =   [0.0000000000000000       -0.0000000000000000       -1.2345678981175834E-006]';
p2 =   [0.0000000000000000        0.0000000000000000        1.0000000000000000]';

for isurf = 1:2
    S(isurf).xroot = readCoeffs2( sprintf( 'c%d.cheb', isurf ) );
    S(isurf).xreg = readCoeffs2( sprintf( 'x%d.bern', isurf ) );
    S(isurf).pn = readCoeffs2( sprintf( 'pn%d.bern', isurf ) );
    S(isurf).sep = reshape( S(isurf).pn, [], 3 );
    S(isurf).sep = importdata( sprintf('sep%d.dat',isurf) );
%     S(isurf).sep = S(isurf).sep ./ repmat( sqrt( sum( S(isurf).sep.^2, 2 ) ), 1, 3 );
end

[ p1b, stat1 ] = separating_plane( S(1).sep, S(2).sep )
[ p2b, stat2 ] = separating_plane( S(1).sep, -S(2).sep )

%%
figure;
hold on

for isurf = 1:2
    if length( S(isurf).sep ) == 3
        quiver3( 0, 0, 0, S(isurf).sep(1), S(isurf).sep(1), S(isurf).sep(1), 0, ...
            'color', 0.5*cls(isurf,:) );
    else
%     bpn = S(isurf).pn;
%     bpn = bpn ./ repmat( sqrt( sum( bpn.^2, 3 ) ), [1,1,3] );
%     surf( bpn(:,:,1), bpn(:,:,2), bpn(:,:,3), ...
%         'facecolor', cls(isurf,:), 'edgecolor', 'k', 'specularstrength', 0 );
    bpn = [ [0,0,0]; S(isurf).sep ];
    [ ep, fp ] = plot_spherical_pyramid( ...
    bpn, convhulln(bpn), 1, 1.0 );
    set( fp, 'facecolor', cls(isurf,:) );
    set( ep, 'color', 'k', 'linewidth', 1 );
    
    sep = importdata( sprintf('sep%d.dat',isurf) );
    size(sep)
    
%     sep = sep ./ repmat( sqrt( sum( sep.^2, 2 ) ), 1, 3 );
    %     plot3( sep(:,1), sep(:,2), sep(:,3), '.', 'color', 0.7*cls(isurf,:) );
%     sep = reshape( sep, [size(S(isurf).pn,1), size(S(isurf).pn,2), 3] );
%     surf( sep(:,:,1), sep(:,:,2), sep(:,:,3), ...
%         'facecolor', 'none', 'edgecolor', 0.5*cls(isurf,:) );
    end
end

[ q, p ] = plot_separating_plane( p1 );

[ q, p ] = plot_separating_plane( p2 );
set( p, 'facecolor', 'g' );

axis image vis3d
view(3)
camlight(30,30);


%%
figure;
hold on

for isurf = 1:2
    creg = chgvar2( S(isurf).xroot, reshape( uvb(isurf,:), 2, 2 ) );
    
    d = cheb_diff2( creg );
    [M,N,~] = size(creg);
    [P,Q] = deal(2*M-1, 2*N-1);
    
    suv = zeros(P,Q,3,2);
    for ivar = 1:2
        for idim = 1:3
            suv(:,:,idim,ivar) = ifcht2d( d(:,:,idim,ivar), P, Q );
        end
    end
    
    sn = zeros(P,Q,3);
    for idim = 1:3
        sn(:,:,idim) = suv(:,:,1+mod(idim,3),1).*suv(:,:,1+mod(idim+1,3),2) - ...
            suv(:,:,1+mod(idim+1,3),1).*suv(:,:,1+mod(idim,3),2) ;
    end

    pn = cat(3, fcht2d(sn(:,:,1)), fcht2d(sn(:,:,2)), fcht2d(sn(:,:,3)) );
    
    si = surf_chebyshev( pn, 1 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
    
    bpn = S(isurf).pn;
    surf( bpn(:,:,1), bpn(:,:,2), bpn(:,:,3), ...
        'facecolor', 'none', 'edgecolor', 0.3*cls(isurf,:) );
    
%     si = surf_chebyshev( creg, 1 );
%     set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
%     
%     bpn = S(isurf).xreg;
%     surf( bpn(:,:,1), bpn(:,:,2), bpn(:,:,3), ...
%         'facecolor', 'none', 'edgecolor', 0.3*cls(isurf,:) );
end

axis image vis3d
view(3)
camlight(30,30);