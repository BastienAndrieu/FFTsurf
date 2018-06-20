clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');


% bc = readCoeffs1( 'root_c_bezier.bern' );


uvbox = [
    0.0000000000000000        1.0000000000000000       -1.0000000000000000        0.0000000000000000     
   0.0000000000000000        1.0000000000000000       -1.0000000000000000        0.0000000000000000];

uvxyz = [
    0.48869071304100015        0.0000000000000000       0.48869071304100004      -0.48869071304100015      -0.42458184934584542        0.0000000000000000       0.27344156094263999     
  0.70322642030743310        0.0000000000000000       0.70322642030743321      -0.70322642030743332      -0.61097370192492040       -2.7755575615628914E-017  0.36135094890504560     
  0.85827630025880253        0.0000000000000000       0.85827630025880242      -0.85827630025880264      -0.74568337209841640        0.0000000000000000       0.44998044005611582     
  0.98360947715858515       -1.6390522841414793E-002  0.96721895431717009       -1.0000000000000000      -0.85457472323820194        1.4240333024641219E-002  0.42633174599902357];


xyz = zeros( size(uvxyz,1), 3, 2 );
for isurf = 1:2
    surfa(isurf).x = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    duv = cheb_diff2( surfa(isurf).x );
    
    reg(isurf).x = chgvar2( surfa(isurf).x, reshape( uvbox(isurf,:), 2, 2 ) );
    reg(isurf).b = readCoeffs2( sprintf('region_%d.bern', isurf) );
    
    surfa(isurf).xu = duv(:,:,:,1);
    surfa(isurf).xv = duv(:,:,:,2);
    
    uv = uvxyz(:,2*isurf + [-1,0]);
    xyz(:,:,isurf) = ICT2unstr( surfa(isurf).x, uv );
    xu = ICT2unstr( surfa(isurf).xu, uv );
    xv = ICT2unstr( surfa(isurf).xv, uv );
    
end


cls = colorcet( 'I2', 'N', 2 );


figure;
hold on

for isurf = 1:2
    %     surf( reg(isurf).b(:,:,1), reg(isurf).b(:,:,2), reg(isurf).b(:,:,3), ...
    %         'facecolor', 'none', 'edgecolor', 0.5*cls(isurf,:) );
    si = surf_chebyshev( surfa(isurf).x, 1 );
    set( si, 'facecolor', cls(isurf,:), 'facealpha', 0.2 );
    
    si = surf_chebyshev( reg(isurf).x, 1 );
    set( si, 'facecolor', cls(isurf,:) );
    
    plot3( xyz(:,1,isurf), xyz(:,2,isurf), xyz(:,3,isurf), ...
        'o', 'color', 0.25*cls(isurf,:), 'markersize', 5 );
end

plot3( uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), 'r*', 'markersize', 5 );

% plot3( bc(:,1), bc(:,2), bc(:,3), 'k-' );

axis image vis3d
view(3)
% view(-130,28)
camproj('persp');
camlight(30,30);




