clc; clear; close all



v = importdata('../trace_intersection_polyline/sep1.dat');
w = importdata('../trace_intersection_polyline/sep2.dat');

% vu = v ./ repmat( sqrt(sum(v.^2,2)), 1, 3 );
% wu = w ./ repmat( sqrt(sum(w.^2,2)), 1, 3 );

'start'
[ q, stat ] = separating_plane(  v, w )


figure; hold on

plot3( v(:,1), v(:,2), v(:,3), 'b.');
plot3( w(:,1), w(:,2), w(:,3), 'r.');

plot3( 0, 0, 0, 'kx');
if stat == 0
    plot_separating_plane( q );
end



p = [-4.2993514424773507E-007  4.9845666108194632E-006  -3.5315297183968447E-006]';
[ qui, pla ] = plot_separating_plane( p );
set( pla, 'facecolor', 'c' );

axis image vis3d
view(3)