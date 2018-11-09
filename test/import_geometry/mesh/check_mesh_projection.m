clc; clear; close all

tri = importdata('tri.dat');
X = [];

for i = 1:2
    figure;
    hold on
    
    switch i
        case 1
            xyz = importdata('xyz.dat');
        case 2
            xyz = importdata('pos_001.dat');
    end
    X = cat(3, X, xyz);
    
    trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
        'facecolor', 'w', 'specularstrength', 0);
    
    axis image %vis3d
    view(3)
    camproj('persp');
    if 1
        [azl,ell] = deal( 120, 60 );
        [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
        light( 'style', 'infinite', 'position', [xl,yl,zl] );
        light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
    end
end


d = sum((X(:,:,1) - X(:,:,2)).^2, 2);
d = 0.5*log10(d);

figure;
trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), d, ...
    'facecolor', 'interp', 'edgecolor', 'none');
axis image %vis3d
view(3)
camproj('persp');
colormap(parula(17));
caxis([-16 0]);
colorbar


