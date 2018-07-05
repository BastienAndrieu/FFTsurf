clc; 
clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');

alpha = 1;

n = 30;

r = [];
for test = 1:100

m = 5000;
k = zeros(m,2);
for i = 1:m
    A = rand(n,n);
    A = 2*A - 1;
%     A = 3*A - 2;
    A = alpha*A;
%     A(:,2) = A(:,1) + (rand/rand) * rand(n,1);    
    
    [Q,R,P] = qr(A);
    
    k(i,:) = [ abs( R(1,1) / R(n,n) ), cond(R) ];
end

% linear regression
ab = polyfit( k(:,1), k(:,2), 1 );
% fprintf('%f\n',ab(1))

r = [ r; ab(1) ];
end
figure; hist(r);
mean(r)
return
% cd(1)
% ef(1)

% sse = sum( (M*ab - k(:,2)).^2 );

% sqrt( 1/m * sum( ( (k(:,2) - polyval(ab,k(:,1)))./k(:,2) ).^2 ) )
% sqrt( 1/m * sum( ( (k(:,2) - ab(1)*k(:,1))./k(:,2) ).^2 ) )

% e = max( abs( k(:,2) - ab(1)*k(:,2) ) );
lo = min( k(:,2)./k(:,1) );
hi = max( k(:,2)./k(:,1) );

mhilo = max( hi, 1/lo );
% lo = 1/mhilo;
% hi = mhilo;
lo = 0.5;
hi = 2;



figure( 'units','pixels','position',[500 200 600 500] );
hold on
plot( k(:,1), k(:,2), '.' )

x = [min(k(:,1)), max(k(:,1))];
% rl1 = plot( x, ab(1)*x, '-', 'linewidth', 2 );

plot( x, lo*x, 'k--', 'linewidth', 1 );
plot( x, hi*x, 'k--', 'linewidth', 1 );

% rl2 = plot( x, x, '--', 'linewidth', 2 );
set( gca, 'yscale', 'log', 'xscale', 'log' );
grid on

% ylabel('cond(A)')
ylabel('$\kappa_2(A)$', ...
    'interpreter', 'latex', 'fontsize', 14, 'rotation', 0, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','right' );
%xlabel('n/2 R_{1,1} / R_{n,n}');
xlabel('$\frac{n}{2} \left| \frac{ R_{11} }{ R_{nn} } \right|$', ...
    'interpreter', 'latex', 'fontsize', 16, ...
    'VerticalAlignment','top' );
title(sprintf('$n = %d$',n), ...
    'interpreter', 'latex', 'fontsize', 14 );
set( gca, 'units', 'normalized', 'position', [0.2 0.15 0.75 0.75 ] )


% legend( rl1, sprintf('$y = %2.2f x$',ab(1)), ...
%     'location', 'northwest', ...
%     'interpreter', 'latex', 'fontsize', 8 );