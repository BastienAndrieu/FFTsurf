clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');


uv0 = [
    0.18416554789070413      -0.98458276729930339       0.96069455607300769      -0.90793532340117189      -0.18416554789071068      -0.98458276729930483
];

%%
for isurf = 1:3
    S(isurf).x = readCoeffs2(sprintf('debug3si_surf%d.cheb',isurf));
end


fid = fopen('debug3si.dat','r');
lowerb = str2num(fgetl(fid));
upperb = str2num(fgetl(fid));
for icurv = 1:2
    curv(icurv).numsurf = str2num(fgetl(fid));
    curv(icurv).np      = str2num(fgetl(fid));
    curv(icurv).uv      = zeros(curv(icurv).np,4);
    for i = 1:curv(icurv).np
        curv(icurv).uv(i,:) = str2num(fgetl(fid));
    end
end
fclose(fid);


%% plot3d
cl = colorcet( 'I2', 'N', 3 );
cl = CC( cl, 0.0, 0.8, 1.5 );


figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

b = zeros(1,3);
for isurf = 1:3
    uvbox = [lowerb(2*isurf-1:2*isurf) ; upperb(2*isurf-1:2*isurf)];
    c = chgvar2(S(isurf).x, uvbox);
%     c = S(isurf).x;
    
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    
    x = ICT2unstr(S(isurf).x, uv0(2*isurf-[1,0]));
    b = b + x/3;
    plot3(x(1), x(2), x(3), 'k.');
end

for icurv = 1:2
    isurf = curv(icurv).numsurf;
    x = ICT2unstr(S(1).x, curv(icurv).uv(:,2*isurf-1:2*isurf));
    plot3( x(:,1), x(:,2), x(:,3), 'k' );
end

% axis image vis3d
axis(b([1,1,2,2,3,3]) + 1e-2*repmat([-1,1],1,3));
daspect([1,1,1])
axis vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


%% plot2d
clr = ['r','b'];




figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
for isurf = 1:3
    subplot(1,3,isurf);
    hold on
    
    rectangle('position',...
        [lowerb(2*isurf-1:2*isurf), ...
        upperb(2*isurf-1:2*isurf) - lowerb(2*isurf-1:2*isurf)]);
    
    if isurf == 1
        lcurv = [1,2];
    else
        lcurv = isurf-1;
    end
    
    for icurv = lcurv
        if isurf == 1
            jsurf = curv(icurv).numsurf;
        else
            jsurf = 1 + mod(curv(icurv).numsurf,2);
        end
        plot(curv(icurv).uv(:,2*jsurf-1), curv(icurv).uv(:,2*jsurf), ...
            '-', 'color', clr(icurv));
        
%         for j = 1:curv(icurv).np
%             text(curv(icurv).uv(j,2*jsurf-1), curv(icurv).uv(j,2*jsurf), ...
%                 num2str(j), 'color', clr(icurv), 'fontsize', 7);
%         end
    end
    
    plot(uv0(2*isurf-1), uv0(2*isurf), 'k*');
    
%     v = u(:,2*isurf-1:2*isurf);
%     for i = 1:size(u,1)-1
%         quiver(v(i,1), v(i,2), v(i+1,1) - v(i,1), v(i+1,2) - v(i,2), 0, 'k');
%     end
    
    axis image
end


%%
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

uv = uv0;
uv = reshape(uv, 2, 3)';
% [ stat, uv, xyz, res, uvpath ] = newton_3_surfaces( S, lowerb, upperb, uv );
[ stat, uv, xyz, res, uvpath ] = newton_3_surfaces_simult_pt_inv( S, lowerb, upperb, uv );
% [ stat, uv, xyz, res, uvpath ] = newton_3_surfaces_min_distance( S, lowerb, upperb, uv );


if 1
    figure;
    plot( res, '.-' );
    set(gca,'yscale','log');
    grid on
    legend({...
        '$|x_1 - x_2|$'; ...
        '$|x_1 - x_3|$'; ...
        '$\approx |u - u_*|$'; ...
        '$\epsilon_{x} \, \mathrm{cond}(\mathbf{J})$'}, ...
        'interpreter', 'latex', 'fontsize', 12, ...
        'location', 'northeast')%outside' )
end




figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
for isurf = 1:3
    subplot(1,3,isurf);
    hold on
    
    rectangle('position',...
        [lowerb(2*isurf-1:2*isurf), ...
        upperb(2*isurf-1:2*isurf) - lowerb(2*isurf-1:2*isurf)]);
    
    plot(uvpath(2*isurf-1,:), uvpath(2*isurf,:), '.-');
    plot(uv(isurf,1), uv(isurf,2), '*');
    
    axis image
end



