addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

%%
tol = 1e-14;

%%
for j = 1:19
    for i = 1:2
        C = readCoeffs2(sprintf('C%d_test%2.2d.txt',i,j));

        [M,N] = cheb_degree_reduction2( C, tol );
        [M,N] = deal(max(3,M),max(3,N));
        C = C(1:M,1:N,:);
        
        D = cheb_diff2( C );
        
        [P,Q] = deal(2*M-1, 2*N-1);
        
        suv = zeros(P,Q,3,2);
        for ivar = 1:2
            for idim = 1:3
                suv(:,:,idim,ivar) = ifcht2d( D(:,:,idim,ivar), P, Q );
            end
        end
        
        sn = zeros(P,Q,3);
        for idim = 1:3
            sn(:,:,idim) = suv(:,:,1+mod(idim,3),1).*suv(:,:,1+mod(idim+1,3),2) - ...
                suv(:,:,1+mod(idim+1,3),1).*suv(:,:,1+mod(idim,3),2) ;
        end
        
        G = cat(3, fcht2d(sn(:,:,1)), fcht2d(sn(:,:,2)), fcht2d(sn(:,:,3)) );
        [P,Q] = cheb_degree_reduction2( G, tol );
        [P,Q] = deal(max(3,P),max(3,Q));
        G = G(1:P,1:Q,:);
        
        writeCoeffs2(G, sprintf('PN%d_test%2.2d.txt',i,j));
    end
end
        
       