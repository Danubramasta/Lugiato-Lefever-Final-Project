function  [eig1,eig2] = dispersion_relation_LLE(X0,Y0,eta,gamma,a,k)

eig1 = - 1 + sqrt(-3*X0^4*eta^2 - 6*X0^2*Y0^2*eta^2 + 4*X0^2*a*eta*k.^2 - 3*Y0^4*eta^2 + 4*Y0^2*a*eta*k.^2 - a^2*k.^4 + 4*X0^2*eta^2*gamma + 4*Y0^2*eta^2*gamma - 2*a*eta*gamma*k.^2 - eta^2*gamma^2);
eig2 = - 1 - sqrt(-3*X0^4*eta^2 - 6*X0^2*Y0^2*eta^2 + 4*X0^2*a*eta*k.^2 - 3*Y0^4*eta^2 + 4*Y0^2*a*eta*k.^2 - a^2*k.^4 + 4*X0^2*eta^2*gamma + 4*Y0^2*eta^2*gamma - 2*a*eta*gamma*k.^2 - eta^2*gamma^2);