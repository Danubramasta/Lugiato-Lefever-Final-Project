clc; clear all; close all;
tic
theta = 0.5;
delta = 1e-3;
Max_iter = 1000;

L = 10;

N = 200;
dx = L/N;

x = [0:dx:(L-dx)]';
D2 = laplacian_1D(N)/dx^2;

E0 = ones(N,1);

eta = 1; a = 1;
gamma = 1;
ttl = ['eta=' num2str(eta) '_a=' num2str(a) '_gamma=' num2str(gamma) '_L=' num2str(L) '_N=' num2str(N) '.mat']%file nya harus .mat tipenya + variabel" diatas

X0 = [real(E0);imag(E0)];

Ei = [0.01 0.011];

[Norm,Ei,lamb,temp_eig,jacobian] = LLE_continuation_fourdif(X0,Ei,eta,gamma,a,N,D2,theta,delta,Max_iter)

save(ttl)
toc