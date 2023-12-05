clc; clear all; close all;
tic
theta = 0.5;
delta = 1e-3;
Max_iter = 1;

L = 50; %gaperlu phi karena phi untuk fungsi yang periodik sedangkan sekarang fungsi nya datar
N = 200;
dx = L/N;

%kebalik woy!!!
m = 2;
[x, D2] = fourdif(N,m); %dengan N banyaknya titik dan m turunan ke berapa
D2 = D2*(2*pi/L)^m; %dinamakan space scalling
x = [0:dx:(L-dx)]';

E0 = ones(N,1);

eta = 1; a = 1;
gamma = 1;
Ei = 1.00;

X0 = [real(E0);imag(E0)]; %tebakan awal ini digunakan untuk fsolve

options = optimset('Display','off','algorithm','levenberg-marquardt');   % Option to display output
[X,fval,exitflag] = fsolve(@(X)LLE_fourdif(X,Ei,eta,gamma,a,N,D2),X0,options);
E = X(1:N)+1i*X(N+1:2*N);
J = jacobian_LLE_fourdif(X,eta,gamma,a,N,D2);
temp_eig = eig(J(1:2*N,1:2*N)); %nilai eigen
lamb = max(real(temp_eig));

% Ei = [0.01 0.011];
% LLE_continuation_fourdif(X0,Ei,eta,gamma,a,N,D2,theta,delta,Max_iter)