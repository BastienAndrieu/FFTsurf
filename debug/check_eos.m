clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

eos_from_faces = 1;
eos_from_edges = 1;
eos_from_verts = 1;

lbshow = 0;
lbscale = 5e-2;
lboffset = 1e-4;
lbuv = -0.9*[1,1];

files_s = dir(fullfile('/d/bandrieu/GitHub/FFTsurf/debug/', ...
    'eos_surf_c_*.cheb'));
files_s = {files_s.name}';
ns = numel(files_s);


files_e = dir(fullfile('/d/bandrieu/GitHub/FFTsurf/debug/', ...
    'eos_edge_c_*.cheb'));
files_e = {files_e.name}';
ne = numel(files_e);

files_v = dir(fullfile('/d/bandrieu/GitHub/FFTsurf/debug/', ...
    'eos_vert_c_*.cheb'));
files_v = {files_v.name}';
nv = numel(files_v);

% cl = colorcet( 'I2', 'N', ns+ne+nv );
% cl = CC( cl, 0.0, 0.8, 1.5 );
% cl = cl(randperm(size(cl,1)),:);
%cl = importdata('/d/bandrieu/GitHub/These/memoire/figures/code/demo_EoS_brep_palette.dat');
cl = importdata('/d/bandrieu/GitHub/These/memoire/figures/data/BRep/face_color.dat');

cl = repmat(cl,2,1);
% cl = ones(ns+ne+nv,3);

cliso = 0.1*cl;

n = 100;
niso = 0;%[9,2];

%%
if 1
    figure;
    hold on
    
    
    pref = 'eos_c_';
    
    if eos_from_faces
        for i = 1:ns
            %         c = readCoeffs2(sprintf('%s%3.3d.cheb',pref,i));
            c = readCoeffs2(files_s{i});
            [si, lu, lv] = surf_chebyshev(c, 1, n, niso);
            set(si, 'facecolor', cl(i,:));
            set([lu;lv], 'color', cliso(i,:));
            
            if lbshow
                plot_local_basis( c, lbuv, lbscale, lboffset );
            end
        end
    end
    
    if eos_from_edges
        for i = 1:ne
            c = readCoeffs2(files_e{i});
            size(c)
            [si, lu, lv] = surf_chebyshev(c, 1, n, niso);
            set(si, 'facecolor', cl(i+ns,:));
            set([lu;lv], 'color', cliso(i+ns,:));
            
            if lbshow
                plot_local_basis( c, lbuv, lbscale, lboffset );
            end
        end
    end
    
    if eos_from_verts
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
    
    axis image vis3d
    view(130,15);
    camlight(-30,30);
    camproj('persp');
end

%%
if 0
    c = readCoeffs2('eos_edge_c_018.cheb');
    seuil = round(log10(eps('double')));
    c = max(seuil,log10(abs(c)));
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
