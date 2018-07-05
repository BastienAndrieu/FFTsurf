clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
% addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');


uv = zeros(2,2,2);
xyz = zeros(3,2);
p = zeros(3,1);
fid = fopen( 'data.dat', 'r' );
fgetl(fid);%
for ipt = 1:2
    for isurf = 1:2
        uv(:,isurf,ipt) = str2num( fgetl(fid) )';
    end
end
fgetl(fid);%
for ipt = 1:2
    xyz(:,ipt) = str2num( fgetl(fid) )';
end
fgetl(fid);%
p = str2num( fgetl(fid) )';
fclose(fid);


x = importdata( 'xyz_polyline.dat' );
u = importdata( 'uv_polyline.dat' );
tng = importdata( 'tangent_polyline.dat' );
xs = tng(:,1:3);
us = tng(:,4:7);

h = sqrt( sum( (x(2:end,:) - x(1:end-1,:)).^2, 2 ) );
s = [ 0 ; cumsum(h) ];
w = x * p;
w = ( w - w(1) ) / ( w(end) - w(1) );



xs*p

% figure;
% plot( w, s / s(end), '.-' );
% grid on
% xlabel('w^*')
% ylabel('s^*')
% return


%% 3d
cls = colorcet( 'I2', 'N', 2 );

figure;
hold on
for isurf = 1:2
    c = readCoeffs2( sprintf('c%d.cheb', isurf) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
    
%     writeCoeffs2( c, sprintf('c%d.cheb', isurf) );
end

plot3( xyz(1,1), xyz(2,1), xyz(3,1), 'ko' );
plot3( xyz(1,2), xyz(2,2), xyz(3,2), 'k*' );
for i = 1:2
    quiver3( xyz(1,i), xyz(2,i), xyz(3,i), p(1), p(2), p(3), 'r' );
end


plot3( x(:,1), x(:,2), x(:,3), 'y.-' );
quiver3( x(:,1), x(:,2), x(:,3), xs(:,1), xs(:,2), xs(:,3), 'b' );


axlim = [ min( [x;xyz'], [], 1 ), max( [x;xyz'], [], 1 ) ];
ctr = 0.5 * ( axlim(4:6) + axlim(1:3) );
rng = 0.5 * ( axlim(4:6) - axlim(1:3) );
rng = 2*rng;
rng = max( 0.5*max(rng), rng );
axlim = ctr([1,1,2,2,3,3]) + rng([1,1,2,2,3,3]) .* repmat([-1,1],1,3);


% axis image
axis(axlim)
daspect( [1,1,1] );
axis vis3d
view(3)
camproj('persp');
camlight(30,30);


% quiver3( xyz(1,1), xyz(2,1), xyz(3,1), ...
%     xyz(1,2) - xyz(1,1), xyz(2,2) - xyz(2,1), xyz(3,2) - xyz(3,1), 0, 'r' );
% 
% y = [1.2392799338878768      -0.23927993388787686       0.14981998347196923];
% plot3( y(1), y(2), y(3), 'b*' );
% 
% return

%% 2d

figure;

for isurf = 1:2
    subplot( 1, 2, isurf )
    hold on
    
    plot( uv(1,isurf,1), uv(2,isurf,1), 'ko' )
    plot( uv(1,isurf,2), uv(2,isurf,2), 'k*' )

        
    plot( u(:,2*isurf-1), u(:,2*isurf), 'k.-' )
    quiver( u(:,2*isurf-1), u(:,2*isurf), us(:,2*isurf-1), us(:,2*isurf), 'b' )
    
%     axis( repmat([-1,1],1,2) );
%     daspect([1,1,1]);
axis image
end
