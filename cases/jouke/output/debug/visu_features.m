clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

plotpaths = 1;


if plotpaths
    verts = importdata('brep/verts.dat');
    vxyz = verts(:,1:3);
    
    fid = fopen('brep/hyperedges.dat', 'r');
    nhe = str2num(fgetl(fid));
    for i = 1:nhe
        hyperedge(i).ne = str2num(fgetl(fid));
        hyperedge(i).verts = str2num(fgetl(fid));
        hyperedge(i).verts
        hyperedge(i).hyperfaces = str2num(fgetl(fid));
        hyperedge(i).edges = zeros(hyperedge(i).ne,2);
        for j = 1:hyperedge(i).ne
            hyperedge(i).edges(j,:) = str2num(fgetl(fid));
        end
    end
    fclose(fid);
    
    fid = fopen('../../init/mesh/paths.dat','r');
    npths = str2num(fgetl(fid));
    for i = 1:npths
        pth(i).hyperedge = str2num(fgetl(fid));
        pth(i).nv = str2num(fgetl(fid));
        pth(i).verts = zeros(1,pth(i).nv);
        for j = 1:pth(i).nv
            l = str2num(fgetl(fid));
            pth(i).verts(j) = round(l(1));
        end
    end
    fclose(fid);
    
    clp = 0.8 * colorcet( 'R3', 'N', max(2,npths) );
else
    idstyp = importdata('../../init/mesh/idstyp.dat');
    l = find(idstyp(:,2) < 2)';
end

instant1 = 0;
instant2 = 52;

% axlim = [-0.2,1.4,-0.6,0.6,-0.6,0.6];
axlim = [0.4,1.1,repmat(0.4*[-1,1],1,2)];

figure('units','normalized', 'position', [.1 .1 .8 .8]);
hold on

for instant = instant1:instant2
    clf;
    hold on
    
    xyz = importdata(sprintf('../pos_%3.3d.dat',instant));
    if plotpaths
        for i = 1:npths
            l = pth(i).verts;
            x = xyz(l,:);
            %         v = hyperedge(pth(i).hyperedge).verts;
            %         x = [vxyz(v(1),:); x; vxyz(v(2),:)];
            plot3(x(:,1), x(:,2), x(:,3), '.-', 'color', clp(i,:), 'linewidth', 1);
        end
    else
        plot3(xyz(l,1), xyz(l,2), xyz(l,3), '.');
    end
    
    axis(axlim);
    axis vis3d
    daspect([1,1,1]);
    view(3)
    camproj('persp');
    
    pause(0.1);
end




















