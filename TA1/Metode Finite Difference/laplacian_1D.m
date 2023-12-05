function D = laplacian_1D(N);

L = ones(1,N-1);
D = -2*eye(N)+diag(L,-1)+diag(L,1);
D(1,2) = 2;
D(N,N-1) = 2;