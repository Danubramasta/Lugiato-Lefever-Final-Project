function F = LLE_fourdif(X,Ei,eta,gamma,a,N,D2)
E = X(1:N)+1i*X(N+1:2*N); %jadi dia vektor ke bawah dari 1:N dia itu bagian realnya,dari N+1:2N itu bagian imajiner
%Exx = ifft((1i*k).^2.*fft(E)); %turunan kedua
Exx = D2*E; %jadiin aray

G = -E+Ei+1i*eta*E.*(abs(E).^2-gamma)+1i*a*Exx;
F = [real(G);imag(G)];
