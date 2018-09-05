clc; clear; close all

figure;
hold on

fid = fopen('debug_curve_segments.dat','r');
while ~feof(fid);
    n = str2num(fgetl(fid));
    x = zeros(n,3);
    for i = 1:n
        x(i,:) = str2num(fgetl(fid));
    end
    plot3( x(:,1), x(:,2), x(:,3), '.-' );
end
fclose(fid);

axis image vis3d
view(3)
camproj('persp');
