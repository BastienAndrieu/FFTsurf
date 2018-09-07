clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/INRIA/Matlab/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf_janvier2016/ImportCAD/');
addpath('/stck8/bandrieu/Bureau/CYPRES/Regression3D/GPeyre/toolbox_graph/');

nfaces = 127;%
cl = colorcet( 'I2', 'N', nfaces );
cl = cl(randperm(size(cl,1)),:);
cl = CC( cl, 0.0, 0.8, 1.5 );

figure('units','normalized','position',[.15 .15 .7 .7]);

hold on

nt = 0;
np = 0;

facets = [];

t = [];
x = [];
ref_t = [];
for iface = 1:nfaces
    tri = importdata(sprintf('tri_%3.3d.dat',iface));
    uv = importdata(sprintf('uv_%3.3d.dat',iface));
%     isurf = importdata(sprintf('numsurf_%3.3d.dat',iface));
%     c = readCoeffs2(sprintf('../propergol/C_%3.3d.cheb',isurf));
    c = readCoeffs2(sprintf('coeffs/c_%3.3d.cheb',iface));
    
    xyz = ICT2unstr(c, uv);
    b = ( xyz(tri(:,1),:) + xyz(tri(:,2),:) + xyz(tri(:,3),:) ) / 3.0;
    
    l = [];%find( b(:,2) < -1.e-4 & b(:,3) > 1.e-4 )';
    l = setdiff(1:size(tri,1), l);
    
    h = trisurf(tri(l,:), xyz(:,1), xyz(:,2), xyz(:,3), ...
        'facecolor', cl(iface,:), 'specularstrength', 0);
    facets = [facets; h];
    
    t = [t ; np+tri];
    x = [x ; xyz];
    ref_t = [ref_t ; iface*ones(size(tri,1),1)];
        
    nt = nt + size(tri,1);
    np = np + size(uv,1);
end



axis image vis3d
view(3)
camproj('persp');
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


%%
[x,t] = fixDuplicateNodes(x, t, 1e-6);
e = compute_edges(t');

inria_mesh( 'fullmesh.mesh', x', t', e, zeros(size(x,1),1), ref_t, zeros(size(e,2),1) );