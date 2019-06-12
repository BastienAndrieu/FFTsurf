clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

lbshow = 0;
lbscale = 5e-2;
lboffset = 1e-4;
lbuv = -0.9*[1,1];

files_s = dir(fullfile('/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/', ...
    'eos_c_*.cheb'));
files_s = {files_s.name}';
ns = numel(files_s);


files_e = dir(fullfile('/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/', ...
    'eos_edge_c_*.cheb'));
files_e = {files_e.name}';
ne = numel(files_e);

files_v = dir(fullfile('/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/', ...
    'eos_vert_c_*.cheb'));
files_v = {files_v.name}';
nv = numel(files_v);

% cl = colorcet( 'I2', 'N', ns+ne+nv );
% cl = CC( cl, 0.0, 0.8, 1.5 );
% cl = cl(randperm(size(cl,1)),:);
cl = importdata('/d/bandrieu/GitHub/These/memoire/figures/code/demo_EoS_brep_palette.dat');
cl = repmat(cl,2,1);
cl = ones(ns+ne+nv,3);

cliso = 0.1*cl;

n = 100;
niso = 0;%[9,2];


if 1
    figure;
    hold on
    
    for j = 2%1:2
        if j == 1
            pref = '../init/coef/c_';
        else
            pref = 'eos_c_';
        end
        
        for i = 1:ns
            c = readCoeffs2(sprintf('%s%3.3d.cheb',pref,i));
            [si, lu, lv] = surf_chebyshev(c, 1, n, niso);
            set(si, 'facecolor', cl(i,:));
            set([lu;lv], 'color', cliso(i,:));
            
            if lbshow
                plot_local_basis( c, lbuv, lbscale, lboffset );
            end
        end
        
        if j == 2
            for i = 1:ne
                c = readCoeffs2(files_e{i});
                [si, lu, lv] = surf_chebyshev(c, 1, n, niso);
                set(si, 'facecolor', cl(i+ns,:));
                set([lu;lv], 'color', cliso(i+ns,:));
                
                if lbshow
                    plot_local_basis( c, lbuv, lbscale, lboffset );
                end
            end
            
            for i = 1:nv
                c = readCoeffs2(files_v{i});
                [si, lu, lv] = surf_chebyshev(c, 1, n, niso);
                set(si, 'facecolor', cl(i+ns+ne,:));
                set([lu;lv], 'color', cliso(i+ns+ne,:));
                
                if lbshow
                    plot_local_basis( c, lbuv, lbscale, lboffset );
                end
            end
        end
    end
    
    axis image vis3d
    view(130,15);
    camlight(-30,30);
    camproj('persp');
end

%%
if 1
%     c = readCoeffs2('eos_c_002.cheb');
    c = readCoeffs2('eos_edge_c_015.cheb');
%     c = readCoeffs2('eos_vert_c_004.cheb');
    c = max(-16,log10(abs(c)));
    cm = parula(17);
    
    figure;
    for i = 1:3
        subplot(1,3,i);
        contourf(c(:,:,i));
        axis equal
        colormap(cm);
        caxis([-16,0]);
    end
end

return
%%
addpath('/stck/bandrieu/Bureau/CYPRES/Propagation/');

figure;
hold on
pref = '../init/coef/c_';
for i = 1:ns
    c = readCoeffs2(sprintf('%s%3.3d.cheb',pref,i));
    
    x = cat(3, ifcht2d(c(:,:,1)), ifcht2d(c(:,:,2)), ifcht2d(c(:,:,3)));
%     s = 0.15*(1.0 + 0.15*cos(5.0*(x(:,:,1)+x(:,:,2)+x(:,:,3))));
    s = 0.15*(1.0 + 0.08*cos(5.0*(x(:,:,1)+x(:,:,2)+x(:,:,3))));
    
    q = EoS_chebyshev_nofrac( x, s, 1, 1 );
    for k = 1:3
        x(:,:,k) = x(:,:,k) + s.*q(:,:,k);
    end
    
    S(i).c = cat(3, fcht2d(x(:,:,1)), fcht2d(x(:,:,2)), fcht2d(x(:,:,3)));
    
    [si, lu, lv] = surf_chebyshev(x, 0, n, niso);
    set(si, 'facecolor', cl(i,:));
    set([lu;lv], 'color', cliso(i,:));
    
    if lbshow
        plot_local_basis( c, lbuv, lbscale, lboffset );
    end
end

axis image vis3d
view(130,15);
camlight(-30,30);
camproj('persp');


c = S(5).c;
c = max(-16,log10(abs(c)));
cm = parula(17);

figure;
for i = 1:3
    subplot(1,3,i);
    contourf(c(:,:,i));
    axis equal
    colormap(cm);
    caxis([-16,0]);
end


