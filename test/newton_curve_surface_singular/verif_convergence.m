clc; clear; %close all

EPS = 2.0*eps('double');
y = importdata( 'residu.dat' );

n = min(200,size(y,1));


order = 1;
mu = order / ( order + 1 );
mrg = 0.25;

% q = log( sum( (y(3:n,4:6)   - y(2:n-1,4:6)).^2, 2) ) ./ ...
%     log( sum( (y(2:n-1,4:6) - y(1:n-2,4:6)).^2, 2) );
q = log( y(3:n,2) ) ./ log( y(2:n-1,2) );

if 0
    figure('name','order of convergence q')
    hold on
    plot( q, 'o-' );
    grid on
    legend( {'$\frac{\log|x_{k+1} - x_k|}{\log|x_k - x_{k-1}|}$'}, ...
        'interpreter', 'latex', 'fontsize', 16 )
end

if mean(q) > 1.5 % *****
    q = 2; % Q-superlinear convergence
else
    q = 1; % linear convergence
end




% sol = importdata('solution.dat');
sol = [ -8.0e-3, 1.0, 0.0 ];
e = sqrt( sum( (y(1:n,4:6) - repmat(sol,n,1)).^2, 2 ) );

if 0
    figure('name','rate of convergence mu')
    hold on
    plot( e.^(1./(1:n)'), '+-' );
    plot( y(3:n,2) ./ y(2:n-1,2), 'x-' );
    grid on
    
    legend( {'$|e_k|^{1/k}$' ; '$\frac{|x_{k+1} - x_k|}{|x_k - x_{k-1}|}$'}, ...
        'interpreter', 'latex', 'fontsize', 14 )
end

normd = y(2:end,2);

figure;
hold on
plot( 0:n-1, y(1:n,1), '.-' );
plot( (1:n-1)-0.0, y(2:n,2), '+-' );
% plot( (1:n-1)-0.5, y(2:n,3)*EPS, '--' );
plot( 0:n-1, y(1:n,3)*EPS, '--' );

% plot( 1:length(normd), normd * ( mu / (1 - mu) ), 'k-' );
% plot( 1:length(normd)-1, normd(1:end-1).*normd(2:end) ./ ( normd(1:end-1) - normd(2:end) ), 'k-' );
% plot( 0:n-1, y(1:n,7), 'k*--' );
% plot( 0:n-2, sqrt( (y(1:n-1,4) - y(n,4)).^2 + (y(1:n-1,5) - y(n,5)).^2 ), 'ko' );
plot( 0:n-1, e, 'ko' );


set( gca, 'yscale', 'log', 'ytick', 10.^(-16:4:4) );%, 'xtick', 0:2:50 );
ylim([1e-16,1]);


legend( {...
    '$|f(x_k)|$' ; ...
    '$|x_{k+1} - x_k|$' ; ...
    '$\epsilon \; \mathrm{cond}\left(\mathbf{J}(x_k)\right)$'; ... '$\approx |x_k - x_*|$'; ...
    '$ |x_k - x_*|$'}, ...
    'interpreter', 'latex', 'fontsize', 12 )
grid on



tuv = y(:,4:6);
figure;
plot( tuv(:,2), tuv(:,1), '.-' );
axis( repmat( [-1,1], 1, 2)  );
daspect([1,1,1])








return


sol = importdata('solution.dat');
rx = sqrt( sum( (y(1:n,4:6) - repmat(sol,n,1)).^2, 2 ) );

figure;
hold on
if order == 0
    plot( rx(2:end) ./ rx(1:end-1).^2, 'o-' );
    %     ylim( 1.6*[1-mrg,1+mrg] );
else
    plot( rx(2:end) ./ rx(1:end-1), 'o-' );
    line( [1,n-1], mu*[1,1], 'color', 'r', 'linestyle', '--' );
    ylim( mu*[1-mrg,1+mrg] );
end
% set( gca, 'yscale', 'log' );
grid on




return





d = y(2:n,2);
% d(2:end) ./ d(1:end-1)

figure;
hold on
if order == 0
    plot( d(2:end) ./ d(1:end-1).^2, 'o-' );
    %     ylim( 1.6*[1-mrg,1+mrg] );
else
    plot( d(2:end) ./ d(1:end-1), 'o-' );
    line( [1,n-1], mu*[1,1], 'color', 'r', 'linestyle', '--' );
    ylim( mu*[1-mrg,1+mrg] );
end
% set( gca, 'yscale', 'log' );
grid on


return
d = y(1:n,1);

figure;
hold on
plot( d(2:end) ./ d(1:end-1), 'o-' );
line( [1,n-1], mu^2*[1,1], 'color', 'r', 'linestyle', '--' );
% set( gca, 'yscale', 'log' );
grid on
ylim( mu^2*[1-mrg,1+mrg] );
