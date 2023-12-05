function J = jacobian_LLE_fourdif(U,eta,gamma,a,N,D2)
X = U(1:N);
Y = U(N+1:2*N);

%notes: kalo diliat kembali, M11-M22 itu blok matriks sedangkan X&Y nya
%vektor maka dari itu M11-M22 harus dijadikan matriks/matriks diagonal 
M11 = diag(-2*X.*Y*eta-1);
M12 = diag(-X.^2*eta-3*Y.^2*eta+eta*gamma)-D2*a;%karna D2 udah matriks makanya gausah diikut diagonalkan
M21 = diag(3*X.^2*eta+Y.^2*eta-eta*gamma)+D2*a; 
M22 = diag(2*X.*Y*eta-1);
%ingat X Y pointwise lagi karna: X Y itu vektor

J = [M11 M12
    M21 M22];