function [ x, feasible ] = lpp_solve( A, c, x )

addpath('/stck/bandrieu/Bureau/CYPRES/LinearProgramming/');

EPSlp = eps('double');
feasible = 1;

debug = 0;

[n,dim] = size(A);
dim = dim - 1;

% A

% if dim == 2
%     A
%     c
%     x
%     visu_lpp( A, c, x )
% end
% fprintf('-----------lpp_solve-----------\n');
% x

if dim == 1
    [ x, feasible ] = solve_LPP_1( A, c );
    return
end

for i = dim+1:n
    %     x
    %     A(i,:)
    
%     fprintf('i = %d, ax + b = %f\n',i, A(i,1:dim)*x + A(i,dim+1) )
    if A(i,1:dim)*x + A(i,dim+1) > (dim+1)*EPSlp; continue; end
%     fprintf('i = %d, x =',i);
%     for k = 1:dim
%         fprintf('%e ',x(k));
%     end
%     fprintf('\ni = %d, Ai =',i);
%     for k = 1:dim+1
%         fprintf('%f ',A(i,k));
%     end
%     fprintf('\n');
%     fprintf('i = %d, ax + b = %e\n',i, A(i,1:dim)*x + A(i,dim+1) )
    
    [~,l] = max( abs(A(i,1:dim)) );
%     fprintf('l = %d\n',l);
    
    if abs(A(i,l)) < 1e-15
        continue
        warning('solve_LPP_d : contrainte principale nulle');
    end
    
    j = setdiff(1:dim,l);
    
    
    Ai = zeros(i-1,dim);
    inv_Ail = 1.0 / A(i,l);
    for k = 1:i-1
        Ai(k,:) = A(k,[j,dim+1]) - A(k,l) * A(i,[j,dim+1]) * inv_Ail;
    end
    ci = c(j) - c(l) * A(i,j)' * inv_Ail;
%     ci = c(j);

    %     xi = x(j) - x(l) * A(i,j)' * inv_Ail;
    xj = x(j);
    
%     ci, Ai, xj
%     keyboard
    if debug
        %         visu_lpp( A(1:i-1,:), c, x )
        l
        visu_lpp( A(1:i,:), c, x )
        Ai, ci, xj
        visu_lpp( Ai, ci, xj )
        %     error('bla')
    end
    
%     keyboard
    % solve lower-dimensional linear programming subproblem
    [ xj, feasible ] = lpp_solve( Ai, ci, xj );
%     xj
    
    if ~feasible
        return
    end
    
    x(j) = xj;
    x(l) = - ( A(i,dim+1) + A(i,j)*xj ) * inv_Ail;
end




end