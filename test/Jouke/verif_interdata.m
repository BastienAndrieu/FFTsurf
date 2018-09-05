clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

nsurftotal = 112;%126;
listsurf = 1:nsurftotal;%[1:7,28:44];%
nsurf = length(listsurf);

PLOT_TREE = 0;
PLOT_NORMALS = 0;

cl = colorcet( 'I2', 'N', nsurftotal );
cl = cl(randperm(size(cl,1)),:);
cl = CC( cl, 0.0, 0.8, 1.5 );
cli = 'k';%'y';

%% read interdata

points = importdata( 'result/interdata_points.dat' );
np = size(points,1);

fid = fopen( 'result/interdata_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );


aabb = 1e6*repmat([1,-1],1,3);
for ic = 1:nc
    curves(ic).dummy = str2num( fgetl( fid ) );
    curves(ic).smooth = str2num( fgetl( fid ) );
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
    
    if n > 0%curves(ic).dummy
        aabb(1:2:5) = min(aabb(1:2:5), min(curves(ic).xyz,[],1));
        aabb(2:2:6) = max(aabb(2:2:6), max(curves(ic).xyz,[],1));
    end
end
fclose( fid );

scn_diameter = norm(aabb(2:2:6) - aabb(1:2:5));


%% Plot

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for isurf = listsurf%1:nsurf%[1:7,28:44]%1:nsurf
    c = readCoeffs2( sprintf('propergol/C_%3.3d.cheb', isurf) );
    surfaces(isurf).D = cheb_diff2(c);
    du = cheb_diff2(surfaces(isurf).D(:,:,:,1));
    dv = cheb_diff2(surfaces(isurf).D(:,:,:,2));
    surfaces(isurf).D2 = cat(4,du,dv(:,:,:,2));
    
    si = surf_chebyshev( c, 1, 100 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
%     set( si, 'facecolor', [1,1,0.7], 'specularstrength', 0 );
    if isurf > 112
        set( si, 'facecolor', 0.5*[1,1,1], 'facealpha', 0.2 );
    end
    
    if PLOT_TREE
        tree = importdata(sprintf('treessi_%d.dat',isurf));
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
    %if curves(ic).dummy; continue; end
    if curves(ic).smooth; continue; end
    if isempty(curves(ic).xyz); continue; end
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
        '-', 'color', cli, 'linewidth', 1.0, 'markersize', 5 )
    e = curves(ic).isplit(:,1)';
    plot3( points(e,1), points(e,2), points(e,3), '.', 'color', cli, 'markersize', 10 );
end

% plot3( points(:,1), points(:,2), points(:,3), '.', 'color', cli, 'markersize', 10 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


if nsurf > 50
    listsurf = [11];%[1:7,28:44]; %
    nsurf = length(listsurf);
    %     return
end

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
if 0
    for ic = 1:nc
        for is = 1:curves(ic).nsplit - 1
            curves(ic).class(is) = 1;
            for jp = 0:1
                ip = curves(ic).isplit(is+jp,1);
                if np2cs(ip) < 2
                    curves(ic).class(is) = -1;
                end
            end
        end
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

l = find(np2cs > 1);
% plot3( points(l,1), points(l,2), points(l,3), '.', 'color', 'b', 'markersize', 10 );

js = 0;
for ic = 1:nc
    %if curves(ic).dummy; continue; end
    if curves(ic).smooth; continue; end
    if isempty(curves(ic).xyz); continue; end
    
    if curves(ic).nsplit > 0
        for is = 1:curves(ic).nsplit - 1
            if curves(ic).class(is) < 0
                clsg = 'none';%[1,0,0];
            else
                clsg = 'b';%[0,1,0];
            end
            
            js = js + 1;
            l = curves(ic).isplit(is,2):curves(ic).isplit(is+1,2);
            plot3( curves(ic).xyz(l,1), curves(ic).xyz(l,2), curves(ic).xyz(l,3), ...
                '-', 'color', clsg, 'markersize', 5 )
            %                 '-', 'color', clrc(js,:), 'markersize', 5 )
            e = curves(ic).isplit(:,1)';
            plot3( points(e,1), points(e,2), points(e,3), '.', 'color', clsg, 'markersize', 10 );
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



%% uv-space
% return
ratio = 1.5;
row = sqrt(nsurf/ratio);
row = max(1, round(row));
col = ceil(nsurf/row);

clcs = CC(cl,0.0,1.2,0.8);

lq = 0.05;%0.05;

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
fid = fopen('result/interdata_surf2curv.dat', 'r');
iisurf = 0;
for isurf=1:nsurftotal
    ns2c = str2num(fgetl(fid));
    if ~isempty(intersect(isurf,listsurf))
        iisurf = iisurf + 1;
        subplot(row,col,iisurf);
        title(sprintf('Surface #%d',isurf));
        hold on
    end
    for jc = 1:ns2c
        icis = str2num(fgetl(fid));
        if ~isempty(intersect(isurf,listsurf))
            ic = icis(1);
            %if curves(ic).dummy; continue; end
            if isempty(curves(ic).xyz); continue; end
            is = icis(2);
            js = icis(3);
            for i = 1:curves(ic).nsplit-1
                if curves(ic).class(i) < 0; continue; end
                l = curves(ic).isplit(i:i+1,2)';
                plot(curves(ic).uv(l(1):l(2),1,is), curves(ic).uv(l(1):l(2),2,is), ...
                    '-', 'color', clcs(js,:), 'linewidth', 1);
                plot(curves(ic).uv(l,1,is), curves(ic).uv(l,2,is), ...
                    '.', 'color', 'k', 'markersize', 7);
                %                 for k = 1:2
                %                     text(curves(ic).uv(l(k),1,is), curves(ic).uv(l(k),2,is), ...
                %                         num2str(curves(ic).isplit(i-1+k,1)));
                %                 end
                
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
        end
    end
    if ~isempty(intersect(isurf,listsurf))
        axis(repmat([-1,1],1,2));
        daspect([1,1,1])
    end
end
fclose(fid);


% return


%% make graphs examples
fid = fopen('result/interdata_surf2curv.dat', 'r');
iisurf = 0;
for isurf=1:nsurftotal
    if ~isempty(intersect(isurf,listsurf))
        iisurf = iisurf + 1;
    end
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
        js = icis(3);
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
            if 0
                pairsurf = [0,0];
                pairsurf(is) = isurf;
                pairsurf(1+mod(is,2)) = js;
                pairsurf
                
                if isempty(intersect(listsurf, js))
                    c = readCoeffs2( sprintf('propergol/C_%3.3d.cheb', js) );
                    surfaces(js).D = cheb_diff2(c);
                    du = cheb_diff2(surfaces(js).D(:,:,:,1));
                    dv = cheb_diff2(surfaces(js).D(:,:,:,2));
                    surfaces(js).D2 = cat(4,du,dv(:,:,:,2));
                end
                
                for k = 1:2
                    [~, duv_ds] = diffgeom_intersection_curve(...
                        surfaces(pairsurf), ...
                        reshape(curves(ic).uv(endpt(k),:,:),2,2));
                    duv_ds
                    duv_ds = (-1)^is * duv_ds(:,is);
                    ang(k) = atan2(duv_ds(2), duv_ds(1));
                end
            else
                if is == 1
                    duv_ds = p([end-1,1],:) - p([end,2],:);
                else
                    duv_ds = p([2,end],:) - p([1,end-1],:);
                end
                ang = atan2(duv_ds(:,2), duv_ds(:,1))';
            end
            
            %             if is == 1
            %                 %                 ang = [atan2(p(end-1,2)-p(end,2), p(end-1,1)-p(end,1)), ...
            %                 %                     atan2(p(1,2)-p(2,2), p(1,1)-p(2,1))];
            %                 ang = atan2(p(1,2)-p(end,2), p(1,1)-p(end,1));
            %             else
            %                 %                 ang = [atan2(p(2,2)-p(1,2), p(2,1)-p(1,1)), ...
            %                 %                     atan2(p(end,2)-p(end-1,2), p(end,1)-p(end-1,1))];
            %                 ang = atan2(p(end,2)-p(1,2), p(end,1)-p(1,1));
            %             end
            %             ang = ang*[1,1];
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
    
    if ~isempty(intersect(isurf,listsurf))
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
    end
    
    %     plot(xy(:,1), xy(:,2), '*');
    %     for iarc = 1:narc
    %         quiver(xy(arc_end(iarc,1),1), xy(arc_end(iarc,1),2), ...
    %             xy(arc_end(iarc,2),1) - xy(arc_end(iarc,1),1), ...
    %             xy(arc_end(iarc,2),2) - xy(arc_end(iarc,1),2), 0 );
    %     end
    %     axis image
end
