function [ x, feasible ] = lp_variant( A, x )
%LP_VARIANT Summary of this function goes here
%   Detailed explanation goes here

% no objective function (any vertex of the feasible region is a solution)
% constraints are in the form of : a(i,1:dim)*x + a(i,dim+1) >= 0

[n,dim] = size(A);
dim = dim - 1;

feasible = 1;

if dim == 1
    [ x, feasible ] = lp_variant_1d( A );
    return
end

for i = dim+1:n
    if A(i,1:dim)*x + A(i,dim+1) >= 0; 
        % the current solution satisfies the i-th constraint
        continue
    end; 
    
    [~,l] = max( abs(A(i,1:dim)) );
    
    j = setdiff(1:dim,l);
    
    Ai = zeros(i-1,dim);
    inv_Ail = 1.0 / A(i,l);
    for k = 1:i-1
        Ai(k,:) = A(k,[j,dim+1]) - A(k,l) * A(i,[j,dim+1]) * inv_Ail;
    end
    
    % solve lower-dimensional linear programming subproblem
    [ xj, feasible ] = lp_variant( Ai, x(j) );
    
    if ~feasible
        return
    end
    
    x(j) = xj;
    x(l) = - ( A(i,dim+1) + A(i,j)*xj ) * inv_Ail;
end

end

