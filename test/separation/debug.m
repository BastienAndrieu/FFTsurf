clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
addpath('/stck/bandrieu/Bureau/CYPRES/SphericalGeometry/');
format long

cls = linspecerp(2);

for i = 1:2
    s(i).x = importdata( sprintf('sep%d.dat',i) );
end


% iaxe = 2;
% for i = 1:2
%     xyi = s(i).x(:,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]);
%     save( sprintf(...
%         '/stck8/bandrieu/Bureau/CYPRES/Geometry/computational_geometry/xy%d.dat', ...
%         8 + i ), ...
%         'xyi', ...
%         '-ascii', '-double' );
% end



[ vec, stat ] = separate_spherical_boxes( s(1).x, s(2).x )


% [ p1, stat1 ] = separating_plane( s(1).x, s(2).x )
% [ p2, stat2 ] = separating_plane( s(1).x, -s(2).x )



figure;
hold on

for i = 1:2
    bpn = [ [0,0,0]; s(i).x ];
    [ ep, fp ] = plot_spherical_pyramid( ...
    bpn, convhulln(bpn), 1, 1.0 );
    set( fp, 'facecolor', cls(i,:) );
    set( ep, 'color', 'k', 'linewidth', 1 );
end

% if stat1 == 0
%     [ q, p ] = plot_separating_plane( p1 );
% end
% 
% if stat2 == 0
%     [ q, p ] = plot_separating_plane( p2 );
%     set( p, 'facecolor', 'g' );
% end

axis image vis3d
view(3)
camlight(30,30);

