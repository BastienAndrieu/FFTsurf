clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

lbshow = 0;
lbscale = 5e-2;
lboffset = 1e-4;
lbuv = -0.9*[1,1];
nisouv = 7;


%%
if 0
    c = readCoeffs2('c_skeleton.cheb');
    
    figure;
    hold on
    
    for k = 1:2
        if k == 2
            for j = 1:size(c,2)
                if mod(j,2) == 0
                    c(:,j,:) = -c(:,j,:);
                end
            end
        end
        
        [s, lu, lv] = surf_chebyshev(c, 1, 100, nisouv);
        set(s, 'facecolor', 'w');
        set([lu;lv], 'color', 'k');
        
        plot_local_basis( c, lbuv, lbscale, lboffset );
    end
    
    axis image vis3d
    view(3)
    camlight(30,30);
    
    return
end


%%
pth = '/d/bandrieu/GitHub/FFTsurf/debug/';
files = dir(fullfile([pth, 'eos_*.cheb']));
files = {files.name}';
ns = numel(files);

cl = colorcet( 'I2', 'N', ns );
cl = CC( cl, 0.0, 0.8, 1.5 );
cl = cl(randperm(size(cl,1)),:);

figure;
hold on

% c = readCoeffs2('c_skeleton.cheb');
% [s, lu, lv] = surf_chebyshev(c, 1, 100, nisouv);
% set(s, 'facecolor', 'w');
% set([lu;lv], 'color', 'k');
% plot_local_basis( c, lbuv, lbscale, lboffset );

for i = 1:ns
    c = readCoeffs2([pth, files{i}]);
    [si, lu, lv] = surf_chebyshev(c, 1, 100, nisouv);
    set(si, 'facecolor', cl(i,:));%, 'facealpha', 0.5);
    set([lu;lv], 'color', 'k');
    
    plot_local_basis( c, lbuv, lbscale, lboffset );
end

axis image vis3d
view(130,15);
camlight(-30,30);
camproj('persp');
