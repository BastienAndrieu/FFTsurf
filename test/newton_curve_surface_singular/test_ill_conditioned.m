clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');

n = 3;
A = rand(n,n);
A(:,2) = A(:,1) + 1e-0 * rand(n,1);

[Q,R,P] = qr(A);

% Determine rank of A
tol = n * eps('double') * abs(R(1,1));
% r = rank(R,tol);
for r = n:-1:1
    if abs(R(r,r)) > tol
        break
    end
end

fprintf('cond = %d,\trank = %d\n', cond(A), r );


b = rand(n,1);
c = Q'*b;

if r == n  % A has full rank
    y = triSup( R(1:n,1:n), c(1:n,:) );
    x = P*y;
    
    e = A*x - b;
    norm(e)
    eps('double')*cond(A)/sqrt(n)
    %while norm(e) > sqrt(n)*eps('double')
    for i = 1:0
%         c = -Q'*e;
%         correc = triSup( R(1:n,1:n), c(1:n,:) );
        correc = -A \ e;
        x = x + correc;
        e = A*x - b;
        norm(e)
    end
    
    
else % A is rank-deficient
    error('rank deficient')
end



