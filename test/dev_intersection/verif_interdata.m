clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

nsurf = 2;
PLOT_TREE = 1;
PLOT_NORMALS = 0;

cl = colorcet( 'I2', 'N', nsurf );
cl = CC( cl, 0.0, 0.8, 1.5 );
cli = 'k';%'y';

%% read interdata

points = importdata( 'interdataglobal_points.dat' );
np = size(points,1);

fid = fopen( 'interdataglobal_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );


aabb = 1e6*repmat([1,-1],1,3);
for ic = 1:nc
    curves(ic).dummy = str2num( fgetl( fid ) );
    %curves(ic).endpoints = str2num( fgetl( fid ) );
    curves(ic).uvbox = zeros( 2, 4 );
    for isurf = 1:2
        curves(ic).uvbox(isurf,:) = str2num( fgetl( fid ) );
    end
    curves(ic).nsplit = str2num( fgetl( fid ) );
    curves(ic).isplit = zeros( curves(ic).nsplit, 2 );
    for ip = 1:curves(ic).nsplit
        curves(ic).isplit(ip,:) = str2num( fgetl( fid ) );
    end
    curves(ic).class = zeros( curves(ic).nsplit-1, 1 );
    for is = 1:curves(ic).nsplit-1
        curves(ic).class(is) = str2num( fgetl( fid ) );
    end
    n = str2num( fgetl( fid ) );
    curves(ic).np = n;
    curves(ic).uv = zeros( n, 2, 2 );
    curves(ic).xyz = zeros( n, 3 );
    for ip = 1:n
        line = str2num( fgetl( fid ) );
        curves(ic).uv(ip,:,1) = line(1:2);
        curves(ic).uv(ip,:,2) = line(3:4);
        curves(ic).xyz(ip,:) = line(5:7);
    end
    
    aabb(1:2:5) = min(aabb(1:2:5), min(curves(ic).xyz,[],1));
    aabb(2:2:6) = max(aabb(2:2:6), max(curves(ic).xyz,[],1));
end
fclose( fid );

scn_diameter = norm(aabb(2:2:6) - aabb(1:2:5));

%% Plot

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for isurf = 1:nsurf
    c = readCoeffs2( sprintf('surfroot%2.2d_x.cheb', isurf) );
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    
    if PLOT_TREE
        tree = importdata(sprintf('treessi_%2.2d.dat',isurf))
        for i = 1:size(tree,1)
            m = round( max(tree(i,[2,4]) - tree(i,[1,3])) * 50 );
            uvi = [ linspace(tree(i,1),tree(i,2),m)' , tree(i,3)*ones(m,1) ;
                tree(i,2)*ones(m,1), linspace(tree(i,3),tree(i,4),m)' ;
                linspace(tree(i,2),tree(i,1),m)' , tree(i,4)*ones(m,1) ;
                tree(i,1)*ones(m,1), linspace(tree(i,4),tree(i,3),m)' ];
            
            l = ICT2unstr( c, uvi );
            
            plot3( l(:,1), l(:,2), l(:,3), '-', ...
                'color', 0.4*cl(isurf,:), 'linewidth', 0.1 );
        end
    end
    
    if PLOT_NORMALS
        d = cheb_diff2(c);
        x = eval_corners_and_center(c);
        xu = eval_corners_and_center(d(:,:,:,1));
        xv = eval_corners_and_center(d(:,:,:,2));
        xn = zeros(5,3);
        for k = 1:5
            xn(k,:) = cross(xu(k,:), xv(k,:));
            xn(k,:) = 5e-2 * scn_diameter * xn(k,:) / norm(xn(k,:));
        end
        quiver3(x(1:4,1), x(1:4,2), x(1:4,3), xn(1:4,1), xn(1:4,2), xn(1:4,3), 0, 'color', 0.3*cl(isurf,:));
    end
end

for ic = 1:nc
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
        '-', 'color', cli, 'linewidth', 1.0, 'markersize', 5 )
end

plot3( points(:,1), points(:,2), points(:,3), '.', 'color', cli, 'markersize', 10 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );



%% point -> curve segment incidence
p2cs = {};
for ip = 1:np
    p2cs{ip} = [];
end

for ic = 1:nc
    for is = 1:curves(ic).nsplit - 1
        for jp = 0:1
            ip = curves(ic).isplit(is+jp,1);
            p2cs{ip} = [ p2cs{ip}; [ic, is] ];
        end
    end
end

np2cs = zeros(np,1);
for ip = 1:np
    np2cs(ip) = size(p2cs{ip},1);
end

%% overwrite curve segments classification
for ic = 1:nc
    for is = 1:curves(ic).nsplit - 1
        curves(ic).class(is) = 1;
        %         for jp = 0:1
        %             ip = curves(ic).isplit(is+jp,1);
        %             if np2cs(ip) < 2
        %                 curves(ic).class(is) = -1;
        %             end
        %         end
    end
end


%%

ncs = 0;
for ic = 1:nc
    if curves(ic).nsplit > 0
        ncs = ncs + curves(ic).nsplit - 1;
    end
end

clrc = hsv(ncs);
clrc = clrc(randperm(ncs),:);

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

plot3( points(:,1), points(:,2), points(:,3), '.', 'color', cli, 'markersize', 10 );

js = 0;
for ic = 1:nc
    if curves(ic).nsplit > 0
        for is = 1:curves(ic).nsplit - 1
            if curves(ic).class(is) < 0
                clsg = [1,0,0];
            else
                clsg = [0,1,0];
            end
            
            js = js + 1;
            l = curves(ic).isplit(is,2):curves(ic).isplit(is+1,2);
            plot3( curves(ic).xyz(l,1), curves(ic).xyz(l,2), curves(ic).xyz(l,3), ...
                '-', 'color', clsg, 'markersize', 5 )
            %                 '-', 'color', clrc(js,:), 'markersize', 5 )
        end
    else
        plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
            '-', 'color', cli, 'markersize', 5 )
    end
end

axis image vis3d
view(3)
camproj('persp');

%%
if 0
    figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
    hold on
    
    for ic = 1:nc
        plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
            '.-', 'color', cli, 'markersize', 5 )
    end
    
    axis image vis3d
    view(3)
    camproj('persp');
end

%%
ratio = 1.5;
row = sqrt(nsurf/ratio);
row = max(1, round(row));
col = ceil(nsurf/row);

clcs = CC(cl,0.0,1.2,0.8);

lq = 0.1;%0.05;

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
fid = fopen('interdataglobal_surf2curv.dat', 'r');
for isurf=1:nsurf
    subplot(row,col,isurf);
    title(sprintf('Surface #%d',isurf));
    hold on
    ns2c = str2num(fgetl(fid));
    for jc = 1:ns2c
        icis = str2num(fgetl(fid));
        ic = icis(1);
        is = icis(2);
        js = icis(3);
        
        for i = 1:curves(ic).nsplit-1
            if curves(ic).class(i) < 0; continue; end
            l = curves(ic).isplit(i:i+1,2)';
            plot(curves(ic).uv(l(1):l(2),1,is), curves(ic).uv(l(1):l(2),2,is), ...
                '-', 'color', clcs(js,:), 'linewidth', 1);
            plot(curves(ic).uv(l,1,is), curves(ic).uv(l,2,is), ...
                '.', 'color', 'k', 'markersize', 7);
            
            p = reshape(curves(ic).uv(l(1):l(2),:,is),[],2);
            ds = sqrt(sum( (p(2:end,:) - p(1:end-1,:)).^2, 2 ));
            s = [0; cumsum(ds)];
            
            if s(end) > 1.5*lq
                if size(p,1) < 3
                    a = 0.5*sum(p,1);
                    if is == 1
                        q = p(1,:) - p(2,:);
                    else
                        q = p(2,:) - p(1,:);
                    end
                else
                    
                    [~,im] = min(abs(s - 0.5*s(end)));
                    a = p(im,:);
                    if is == 1
                        q = p(im-1,:) - a;
                    else
                        q = p(im+1,:) - a;
                    end
                end
                q = lq * q / norm(q);
                w = 0.3*[q(2), -q(1)];
                h = repmat(a,3,1) + [w;q;-w];
                patch('xdata', h(:,1), 'ydata', h(:,2), ...
                    'edgecolor', 'none', 'facecolor', clcs(js,:));
            end
        end
        %         l = curves(ic).isplit(:,2)';
        %         plot(curves(ic).uv(l,1,is), curves(ic).uv(l,2,is), '.', 'color', 0.7*cl(isurf,:), 'markersize', 10);
        %
        %         plot(curves(ic).uv(:,1,is), curves(ic).uv(:,2,is), '-', 'color', 0.7*cl(isurf,:));
    end
    axis(repmat([-1,1],1,2));
    daspect([1,1,1])
end
fclose(fid);


return


%% make graphs examples
fid = fopen('interdataglobal_surf2curv.dat', 'r');
for isurf=1:nsurf
    
%     figure;
%     hold on
    
    arc_end = [];
    arc_ang = [];
    xy = [];
    ns2c = str2num(fgetl(fid));
    for jc = 1:ns2c
        icis = str2num(fgetl(fid));
        ic = icis(1);
        is = icis(2);
        
        if is == 1
            order = [2,1];
        else
            order = [1,2];
        end
        
        for i = 1:curves(ic).nsplit-1
            endpt = curves(ic).isplit(i:i+1,2)';
            endve = curves(ic).isplit(i:i+1,1)';
            
            p = reshape(curves(ic).uv(endpt(1):endpt(2),:,is),[],2);
            %             plot(p([1,end],1), p([1,end],2), 'o');
            
            xy(endve(1),:) = p(1,1:2);
            xy(endve(2),:) = p(end,1:2);
            
            %             quiver(xy(endve(order(1)),1), xy(endve(order(1)),2), ...
            %                 xy(endve(order(2)),1) - xy(endve(order(1)),1), ...
            %                 xy(endve(order(2)),2) - xy(endve(order(1)),2), 0);
            
            arc_end = [arc_end; endve(order)];
            
            if is == 1
                %                 ang = [atan2(p(end-1,2)-p(end,2), p(end-1,1)-p(end,1)), ...
                %                     atan2(p(1,2)-p(2,2), p(1,1)-p(2,1))];
                ang = atan2(p(1,2)-p(end,2), p(1,1)-p(end,1));
            else
                %                 ang = [atan2(p(2,2)-p(1,2), p(2,1)-p(1,1)), ...
                %                     atan2(p(end,2)-p(end-1,2), p(end,1)-p(end-1,1))];
                ang = atan2(p(end,2)-p(1,2), p(end,1)-p(1,1));
            end
            ang = ang*[1,1];
            arc_ang = [arc_ang; ang];
        end
    end
    l = unique(arc_end(:)');
    
    linv = zeros(1,size(xy,1));
    for i = 1:length(l)
        linv(l(i)) = i;
    end
    xy = xy(l,:);
    arc_end = linv(arc_end);
    
    narc = size(arc_end,1);
    nnod = size(xy,1);
    
    fout = fopen(sprintf('../graph/test%2.2d.dat',isurf),'w');
    fprintf(fout,'%d\n',narc);
    for iarc = 1:narc
        fprintf(fout, '%d %d\n', arc_end(iarc,1), arc_end(iarc,2));
    end
    for iarc = 1:narc
        fprintf(fout, '%f %f\n', arc_ang(iarc,1), arc_ang(iarc,2));
    end
    fprintf(fout,'%d\n',nnod);
    for inod = 1:nnod
        fprintf(fout, '%f %f\n', xy(inod,1), xy(inod,2));
    end
    fclose(fout);
    
%     plot(xy(:,1), xy(:,2), '*');
%     for iarc = 1:narc
%         quiver(xy(arc_end(iarc,1),1), xy(arc_end(iarc,1),2), ...
%             xy(arc_end(iarc,2),1) - xy(arc_end(iarc,1),1), ...
%             xy(arc_end(iarc,2),2) - xy(arc_end(iarc,1),2), 0 );
%     end
%     axis image
end


