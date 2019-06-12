clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Propagation/');

normal_speed = @(x,y,z) 1. + 0.15*cos(5.*(x+y+z));

pth = '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/';
files = dir(fullfile([pth, 'eos_*.cheb']));
files = {files.name}';
ns = numel(files);

for i = 1:ns
    S(i).c = readCoeffs2(files{i});
end

timestep = 1.e-2;
ninstants = 40;

for instant = 1:ninstants
    fprintf('instant %d/%d\n',instant,ninstants);
    for i = 1:ns
        xyz = cat(3, ifcht2d(S(i).c(:,:,1)), ifcht2d(S(i).c(:,:,2)), ifcht2d(S(i).c(:,:,3)));
        speed = normal_speed(xyz(:,:,1), xyz(:,:,2), xyz(:,:,3));
        q  = EoS_chebyshev_nofrac( xyz, speed, timestep, 1 );
        for j = 1:3
            xyz(:,:,j) = xyz(:,:,j) + timestep*speed.*q(:,:,j);
        end
        S(i).c = cat(3, fcht2d(xyz(:,:,1)), fcht2d(xyz(:,:,2)), fcht2d(xyz(:,:,3)));
    end
end


%%
cl = colorcet( 'I2', 'N', ns );
cl = CC( cl, 0.0, 0.8, 1.5 );
cl = cl(randperm(size(cl,1)),:);

figure;
hold on

for i = 1:ns
    [si, lu, lv] = surf_chebyshev(S(i).c, 1, 100);
    set(si, 'facecolor', cl(i,:));
end

axis image vis3d
view(130,15);
camlight(-30,30);
camproj('persp');