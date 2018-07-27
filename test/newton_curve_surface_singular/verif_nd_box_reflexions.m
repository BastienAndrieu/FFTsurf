clc; clear; close all
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

EPS = 10*eps('double');


l = [-0.50000000000000000      -0.12500000000000000       0.16535034245560629];
u = [-0.48332957322768932       -9.4986984464060190E-002  0.15625000000000000];

x = [-0.49166478661384466      -0.10999349223203009       0.16080017122780316];
d = [-6.6720644512181220E-003   7.5035799421662320E-003   4.5590883763657825E-003];

dim = length(x);

% check if x is inside feasible region...
for idim = 1:dim
    if x(idim) < l(idim) - EPS || x(idim) > u(idim) + EPS
        fprintf('!!! start point out of bounds\n');
        break
    end
end


%%
xtmp = x;
dxtmp = d;

reflmax = 1000;
for irefl = 1:reflmax
    lambda = 1.0;
    ib = 0;
    
    for idim = 1:dim
        if ( abs(dxtmp(idim)) < EPS ), continue; end
        if ( dxtmp(idim) < 0.0 )
            lambda_i = min( lambda, (l(idim) - xtmp(idim)) / dxtmp(idim) );
        else
            lambda_i = min( lambda, (u(idim) - xtmp(idim)) / dxtmp(idim) );
        end
        if ( lambda > lambda_i )
            ib = idim;
            lambda = lambda_i;
        end
    end
    
    if ( ib == 0 )
        xnew = xtmp + dxtmp;
        dx = xtmp - x;
        return
    else
        xtmp = xtmp + lambda*dxtmp;
        dxtmp = (1.0 - lambda) * dxtmp;
        dxtmp(ib) = -dxtmp(ib);
    end
end
           
%% plot
boxcolor = 'k';
boxtolcolor = 'k';%'r';
xnewcolor = 'b';

figure;
hold on

h = plotAABB( [l(1), u(1), l(2), u(2), l(3), u(3)] );
set( h, 'facecolor', boxcolor, 'facealpha', 0.1, 'edgecolor', boxcolor );


plot3( x(1), x(2), x(3), 'ko' );

quiver3( x(1), x(2), x(3), d(1), d(2), d(3), 0, 'r' );

% plot3( xnew(1), xnew(2), xnew(3), '*', 'color', xnewcolor );

view(3)
axis equal vis3d
%     camlight(30,30);
