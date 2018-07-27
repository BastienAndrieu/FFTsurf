clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');

cl = colorcet( 'I2', 'N', 2 );


[icurv, ivar, ival] = deal(1,1,1);

uvbox = [
    0.0000000000000000       0.50000000000000000       -1.0000000000000000      -0.50000000000000000     
  -1.0000000000000000      -0.50000000000000000       -1.0000000000000000      -0.50000000000000000
    ];

tbox = [
    -1,1
    ];
ubox = [
    -1.0000000000000000      -0.50000000000000000       -1.0000000000000000      -0.50000000000000000
];

for isurf = 1:2
    c = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    if isurf ~= icurv
        S(isurf).x = chgvar2(c, reshape(ubox, 2, 2));
    else
        S(isurf).x = chgvar2(c, reshape(uvbox(isurf,:), 2, 2));
    end
    
    d = cheb_diff2(S(isurf).x);
    S(isurf).xu = d(:,:,:,1);
    S(isurf).xv = d(:,:,:,2);
    
    if isurf == icurv
        C.x = border_parameterization(S(isurf).x, ivar, ival);
        C.x = chgvar1(C.x, tbox);
        C.xt = cheb_diff1(C.x);
    end
end


n = 200;
t = linspace(-1, 1, n)';

figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on
for isurf = 1:2%1+mod(icurv,2)%
    si = surf_chebyshev(S(isurf).x, 1);
    set(si, 'facecolor', cl(isurf,:), 'specularstrength', 0);
end

c = chebval1(C.x, t);
plot3(c(:,1), c(:,2), c(:,3), 'k-');

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);

uvxyz = [
0.0000000000000000      -0.30494613293843287      -0.10776183717734601       -4.9314229082246912E-002  -4.7846698502724474E-017   1.4877639051280589E-014   1.6028844918025698E-015
   4.5597736279888836E-003 -0.19761684264029328        0.0000000000000000       -4.7024838452958018E-002  -5.7169249364045042E-003  -6.7385879221049133E-002  -2.2893906292872773E-003
  0.10012042262145293      -0.76640502470232408      -0.61467581337170407        0.0000000000000000      -0.12024570621814414       0.28589777038412822       -4.9314229082245260E-002
  ];

plot3(uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), 'k.');


return
writeCoeffs2(S(1+mod(icurv,2)).x,'../newton_curve_surface_singular/surf.cheb');
writeCoeffs1(C.x,'../newton_curve_surface_singular/curv.cheb');

%%
% bc = chebyshev2bezier_1(C.x);
% bs = chebyshev2bezier_2(S(1+mod(icurv,2)).x);
bc = readCoeffs1('debugcsstc_bc.bern');
bs = readCoeffs2('debugcsstc_bs.bern');


[M,N,~] = size(bs);
P = size(bc,1);

for k = 1:2
    for j = 1:2
        s = reshape(bs((j-1)*(M-1)+1,(k-1)*(N-1)+1,:), 1, 3);
        for i = 1:2
            c = bc((i-1)*(P-1)+1,:);
            if norm(s - c) < 1e-6
                t = (-1)^i;
                uv = (-1).^[j,k]
%                 t = tbox(i);
%                 uv = ubox([j,2+k]);
                x = 0.5*(c + s);
            end
        end
    end
end
% uv = [-0.10776183717734601       -4.9314229082246912E-002];
x = [-4.7846698502724474E-017   1.4877639051280589E-014   1.6028844918025698E-015];

s = ICT2unstr(S(1+mod(icurv,2)).x, uv);
su = ICT2unstr(S(1+mod(icurv,2)).xu, uv);
sv = ICT2unstr(S(1+mod(icurv,2)).xv, uv);
n = cross(su, sv)'

n = n / norm(n);

plot3(s(1), s(2), s(3), '*');

quiver3(x(1), x(2), x(3), su(1), su(2), su(3), 'r');
quiver3(x(1), x(2), x(3), sv(1), sv(2), sv(3), 'g');
quiver3(x(1), x(2), x(3), n(1), n(2), n(3), 'b');


seps = reshape(bs,[],3);
seps = seps - repmat(x, size(seps,1), 1);

sepc = bc;
sepc = sepc - repmat(x, size(sepc,1), 1);

dotns = seps*n;
dotnc = sepc*n;

for i = 1:P
    fprintf('%d \t %d\n',i,dotnc(i));
end

for j = 1:N
    for i = 1:M
        fprintf('%d %d \t %d\n',i,j,dotns((j-1)*M+i));
    end
end

% [min(dotns), max(dotns)]
% [min(dotnc), max(dotnc)]

max( [min(dotns), min(dotnc)] ) - min( [max(dotns), max(dotnc)] )

return
%%

ct = chebval1(C.xt, t)';
ct = ct / norm(ct);
[a,b] = completeBase(ct);
rot = [a,b,ct];

EPSxyz = 1e-9;

seps = reshape(bs,[],3);
seps = seps - repmat(x, size(seps,1), 1);
seps(sqrt(sum(seps.^2,2)) < EPSxyz,:) = [];
seps = seps ./ repmat(sqrt(sum(seps.^2,2)), 1, 3);

sepc = bc;
sepc = sepc - repmat(x, size(sepc,1), 1);
sepc(sqrt(sum(sepc.^2,2)) < EPSxyz,:) = [];
sepc = sepc ./ repmat(sqrt(sum(sepc.^2,2)), 1, 3);


quiver3(0, 0, 0, a(1), a(2), a(3), 'r');
quiver3(0, 0, 0, b(1), b(2), b(3), 'g');
quiver3(0, 0, 0, ct(1), ct(2), ct(3), 'b');



sepc = sepc*rot;
seps = seps*rot;

seps(abs(abs(seps(:,3))-1) < EPSxyz,:) = [];
% abs(abs(sepc(:,3))-1)
sepc(abs(abs(sepc(:,3))-1) < EPSxyz,:) = [];


[ vec, stat ] = separate_spherical_boxes(sepc, seps, 3);


figure;
hold on
plot(0,0,'kx');
q = convhulln(seps(:,1:2));

% wedgec = minimal_2d_wedge_2( sepc(:,1:2) )
% wedges = minimal_2d_wedge_2( seps(:,1:2) )
% return
% pw = plot_wedge( wedges );
%     set( pw, 'facecolor', 'c', 'edgecolor', 'none' );
% pw = plot_wedge( wedgec );
%     set( pw, 'facecolor', 'y', 'edgecolor', 'none' );
    
plot(seps(:,1), seps(:,2), '.', 'color', 0.7*cl(1+mod(icurv,2),:));
for i = 1:size(q,1)
%     plot(seps(q(i,:),1), seps(q(i,:),2), 'r-');
end
plot(sepc(:,1), sepc(:,2), '.', 'color', 0.7*cl(icurv,:));
axis( repmat([-1,1],1,2) );
daspect([1,1,1])




% [ vec, stat ] = separating_plane( sepc, seps )

% qs = myconvhulln(seps);
% qc = myconvhulln(sepc);

figure;
hold on

plot3(0,0,0,'kx');

% surf(S.b(:,:,1)-tuvxyz(1,4), S.b(:,:,2)-tuvxyz(1,5), S.b(:,:,3)-tuvxyz(1,6), ...
%     'facecolor', 'none', 'edgecolor', 0.5*cl(1+mod(icurv,2),:));
% plot3(G.b(:,1)-tuvxyz(1,4), G.b(:,2)-tuvxyz(1,5), G.b(:,3)-tuvxyz(1,6), ...
%     '.-', 'color', 0.5*cl(icurv,:));

for isurf = 1:2
    if isurf == icurv
%         q = qc;
        sep = sepc;
    else
%         q = qs;
        sep = seps;
    end
    plot3( sep(:,1), sep(:,2), sep(:,3), ...
        '.', 'color', 0.7*cl(isurf,:));
%     trisurf(q, sep(:,1), sep(:,2), sep(:,3), ...
%         'facecolor', cl(isurf,:), 'edgecolor', 'none', 'facealpha', 1 );
end

if stat == 0
    plot_separating_plane(vec);
end

% axis image vis3d
axis( repmat([-1,1],1,3) );
daspect([1,1,1])
axis vis3d
view(120,30)
camproj('persp');
camlight(30,30);
    

return
%%
EPSuv = 0;
lowerb = reshape(uvbox(:,[1,3]), 4, 1) - EPSuv;
upperb = reshape(uvbox(:,[2,4]), 4, 1) + EPSuv;

uv = 0.5*(uvbox(:,[1,3]) + uvbox(:,[2,4]));
% uv(ivar,icurv) = uvbox(isurf,2*(ivar-1)+ival);


for it = 1:20
    
    for isurf = 1:2
        xyz(:,isurf) = ICT2unstr(S(isurf).x, uv(isurf,:))';
    end
    r = xyz(:,1) - xyz(:,2);
    
    for isurf = 1:2
        jac(:,2*isurf-1) = ICT2unstr(S(isurf).xu, uv(isurf,:))';
        jac(:,2*isurf)   = ICT2unstr(S(isurf).xv, uv(isurf,:))';
    end
    jac(:,3:4) = -jac(:,3:4);
    
    duv = -jac\r;
    res(it,:) = [norm(r),norm(duv),cond(jac)*eps('double')];
    
    %     duv = nd_box_reflexions(uv(:), lowerb, upperb, duv);
    
    
    for isurf = 1:2
        uv(isurf,:) = uv(isurf,:) + duv(2*isurf-1:2*isurf)';
    end
end

for isurf = 1:2
    xyz(:,isurf) = ICT2unstr(S(isurf).x, uv(isurf,:))';
    plot3(xyz(1,isurf), xyz(2,isurf), xyz(3,isurf), '*', 'color', 0.7*cl(isurf,:));
end
xyz = 0.5*sum(xyz,2);
plot3(xyz(1), xyz(2), xyz(3), 'ro');


figure;
plot(res, '.-');
set(gca, 'yscale', 'log');
grid on
legend({'|x_1 - x_2|','|\Delta u|','tol'});