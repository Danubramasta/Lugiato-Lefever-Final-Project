function F = LLE(X,Ei,eta,gamma,a,N,k)
E = X(1:N)+1i*X(N+1:2*N);
Exx = ifft((1i*k).^2.*fft(E)); 

G = -E+Ei+1i*eta*E.*(abs(E).^2-gamma)+1i*a*Exx;
F = [real(G);imag(G)];
