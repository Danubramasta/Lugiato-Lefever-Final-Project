clc; clear all; close all;
L = 3*pi;
% 
 N = 200;
% dx = L/N;
% 
% x = [-L/2:dx:(L/2-dx)]';
% 
kx = [0:N/2-1 -N/2:-1]'*2*pi/L;

E0 = ones(N,1);
Ei = 1.5;%0.01;
eta = 1; a = 1;
gamma = 1;
X0 = [real(E0);imag(E0)];

options = optimset('Display','iter','algorithm','levenberg-marquardt');   % Option to display output
[X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE(X,Ei,eta,gamma,a,N,kx),X0,options);
E = X(1:N)+1i*X(N+1:2*N);

jac = jacobian(1:2*N,1:2*N);
temp_eig = eig(jacobian(1:2*N,1:2*N));
lamb1 = max(real(temp_eig));

% figure(2)
% plot(real(temp_eig),imag(temp_eig),'o')

p = -10:0.1:10;
X0 = X(1); Y0 = X(N+1);
[eig1,eig2] = dispersion_relation_LLE(X0,Y0,eta,gamma,a,p);

figure(1)
subplot(2,1,1)
plot(p,real(eig1),p,real(eig2))
% xlabel('p')
ylabel('Re(\lambda)')
title(['Kestabilan Bagian Real'])

subplot(2,1,2)
plot(p,imag(eig1),p,imag(eig2))
%xlabel('\lambda max')
ylabel('Im(\lambda)')
title(['Kestabilan Bagian Imajiner'])

kc = sqrt(-a*eta*(-2*X0^2 - 2*Y0^2 + gamma))/a;
[eig1,eig2] = dispersion_relation_LLE(X0,Y0,eta,gamma,a,kc);

lamb = eig1;

[lamb1 lamb]
