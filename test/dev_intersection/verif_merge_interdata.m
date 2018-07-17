clc; clear; close all

cl = linspecer(3);
ls = {'-';'-';'--'};

%%
figure;
hold on

fid = fopen('merge_interdata.dat','r');
for i = 1:5
    b = str2num( fgetl(fid) );
    if i <= 3
        rectangle( 'position', [b([1,3]), b([2,4]) - b([1,3])], ...
            'edgecolor', cl(i,:), 'linestyle', ls{i} );
    end
    if i >= 3
        uvbox(i-2,:) = b;
    end
end

for i = 1:2
    n = str2num( fgetl(fid) );
    u = zeros(n,2);
    x = zeros(n,3);
    for j = 1:n
        line = str2num( fgetl(fid) );
        u(j,:) = line(1:2);
        v(j,:) = line(3:4);
        x(j,:) = line(5:7);
    end
    plot( u(:,1), u(:,2), '.-', 'color', cl(i,:) );
    
    curv(i).n = n;
    curv(i).uv1 = u;
    curv(i).uv2 = v;
    curv(i).x = x;
end
fclose( fid );

axis image




%%
nsurf = 3;

cl = colorcet( 'I2', 'N', nsurf );
cl = CC( cl, 0.0, 0.8, 1.5 );

% surfaces and curves
figure;
hold on

for isurf = 1:nsurf
    %     c = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    c = readCoeffs2( sprintf('debugmrg_surf%d.cheb', isurf) );
    S(isurf).x = c;
    c = chgvar2(c, reshape(uvbox(isurf,:),2,2));
    
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
end

for i = 1:2
    x = curv(i).x;
    plot3( x(:,1), x(:,2), x(:,3), 'k' );
end

axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


% curves only
figure;
hold on

for i = 1:2
    x = curv(i).x;
    plot3( x(:,1), x(:,2), x(:,3), 'k' );
end

axis image vis3d
view(3)
camproj('persp');


%% intersect polylines in uv space of surface 1
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/polylines/');
[ isegm, t ] = intersect_polylines( curv(1).uv1, curv(2).uv1, [1,1], [curv(1).n,curv(2).n], [], [] )

uv = zeros(3,2);
uv(1,:) = (1.0 - t(1)) * curv(1).uv1(isegm(1),:) + t(1) * curv(1).uv1(isegm(1)+1,:);
uv(2,:) = (1.0 - t(1)) * curv(1).uv2(isegm(1),:) + t(1) * curv(1).uv2(isegm(1)+1,:);
uv(3,:) = (1.0 - t(2)) * curv(2).uv2(isegm(2),:) + t(2) * curv(2).uv2(isegm(2)+1,:);

uv

%% Newton refinement
for isurf = 1:nsurf
    d = cheb_diff2( S(isurf).x );
    S(isurf).xu = d(:,:,:,1);
    S(isurf).xv = d(:,:,:,2);
end

EPSxyz = 1e-9;
EPSuv  = 1e-12;
itmax  = 20;

xyz   = zeros(3,3);
jac   = zeros(6,6);
duv   = zeros(6,1);
erruv = 0.0;
condJ = 1.0;
res = [];
for it = 1:itmax
    % compute residual
    for isurf = 1:3
        xyz(:,isurf) = ICT2unstr(S(isurf).x, uv(isurf,:))';
    end
    
    r12 = xyz(:,1) - xyz(:,2);
    r13 = xyz(:,1) - xyz(:,3);
    
    res = [ res ; [ norm(r12), norm(r13), erruv, max(EPSuv, eps('double')*condJ) ] ];
    
    % termination criteria
    if erruv < max(EPSuv, eps('double')*condJ)
        if norm(r12) < EPSxyz && norm(r13) < EPSxyz
            xyz = sum(xyz, 2) / 3.0;
            break
        end
    end
    
    % compute Jacobian matrix
    jac(1:3,1) =  ICT2unstr(S(1).xu, uv(1,:))';
    jac(1:3,2) =  ICT2unstr(S(1).xv, uv(1,:))';
    jac(1:3,3) = -ICT2unstr(S(2).xu, uv(2,:))';
    jac(1:3,4) = -ICT2unstr(S(2).xv, uv(2,:))';
    jac(4:6,1) =  jac(1:3,1);
    jac(4:6,2) =  jac(1:3,2);
    jac(4:6,5) = -ICT2unstr(S(3).xu, uv(3,:))';
    jac(4:6,6) = -ICT2unstr(S(3).xv, uv(3,:))';
    
    % solve for Newton step
    [U,W,V] = svd(jac);
    wmax = max(diag(W));
    wmin = min(diag(W));
    condJ = wmax/wmin;
    wmin = 6.0 * eps('double') * wmax;
    for i = 1:6
        if W(i,i) > wmin
            W(i,i) = 1.0 / W(i,i);
        else
            W(i,i) = 0;
        end
    end
    
    duv = - V * W * U' * [r12 ; r13];
    
    % keep iterate in feasible domain
    % ...
    
    % update solution
    for isurf = 1:3
        uv(isurf,:) = uv(isurf,:) + duv(2*isurf-1:2*isurf)';
    end
end

if 1
    figure;
    plot( res, '.-' );
    set(gca,'yscale','log');
    grid on
end


return

%% Intersect elsewhere?
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/SphericalGeometry/');

for isurf = 1:nsurf
    c = chgvar2(S(isurf).x, reshape(uvbox(isurf,:),2,2));
    S(isurf).b = chebyshev2bezier_2(c);
    
    st = ( uv(isurf,:) - uvbox(isurf,[1,3]) ) ./ ( uvbox(isurf,[2,4]) - uvbox(isurf,[1,3]) );
    
    [ bsw, bse, bnw, bne ] = subdivide_Bezier_patch( S(isurf).b, st );
    
    S(isurf).child(1).b = bsw;
    S(isurf).child(2).b = bse;
    S(isurf).child(3).b = bnw;
    S(isurf).child(4).b = bne;
    
    colhsv = repmat( rgb2hsv( cl(isurf,:) ), 4, 1 );
    colhsv(:,1) = mod( colhsv(:,1) + 0.07*(2*rand(4,1)-1), 1.0 );
    colhsv(:,2) = min( max( colhsv(:,2) + 0.2*(2*rand(4,1)-1), -1.0 ), 1.0 );
    colhsv(:,3) = min( max( colhsv(:,3) + 0.3*(2*rand(4,1)-1), -1.0 ), 1.0 );
    cl = [ cl; hsv2rgb( colhsv ) ];
end


for isurf = 1:nsurf
    for ichild = 1:4
        bcp = S(isurf).child(ichild).b;
        for ivar = 1:3
            bcp(:,:,ivar) = bcp(:,:,ivar) - xyz(ivar);
        end
        sep = reshape(bcp, [], 3);
        sep(sqrt(sum(sep.^2,2)) < EPSxyz,:) = [];
        S(isurf).child(ichild).sep = sep ./ repmat(sqrt(sum(sep.^2,2)), 1, 3);
    end
end


n = 100;
t = linspace(0,1,n)';

for kchild = 1:4
    for jchild = 1:4
        for ichild = 1:4
            
            ic = [ichild,jchild,kchild];
            
            [~,stat] = separating_plane(S(1).child(ichild).sep, S(2).child(jchild).sep);
            if stat == 0; continue; end;
            [~,stat] = separating_plane(S(1).child(ichild).sep, S(3).child(kchild).sep);
            if stat == 0; continue; end;
            [~,stat] = separating_plane(S(2).child(jchild).sep, S(3).child(kchild).sep);
            if stat == 0; continue; end;
            
            ic
            
            figure;
            hold on
            
            for isurf = 1:nsurf
                s = eval_Bezier_patch_tpgrid(S(isurf).child(ic(isurf)).b, t, t);
                surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                    'edgecolor','none','facecolor',cl(3+4*(isurf-1)+ichild,:),'diffusestrength',0.7, ...
                    'facelighting','gouraud','specularstrength',0);
                %                 sep = [ [0,0,0] ; S(isurf).child(ic(isurf)).sep ];
                %                 q = convhulln(sep);
                %                 [eps,fps] = plot_spherical_pyramid( sep, q, 1, 1 );
                %                 set( eps, 'color', 'k', 'linewidth', 1.0 );
                %                 set( fps, 'facecolor', cl(3+4*(isurf-1)+ichild,:), 'edgecolor', 'k' );
            end
            
            
            axis image vis3d
            view(3)
            camproj('persp');
            % camlight(30,30);
            [azl,ell] = deal( 120, 60 );
            [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
            light( 'style', 'infinite', 'position', [xl,yl,zl] );
            light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
            
        end
    end
end

