clc; clear all; close all;
tic
theta = 0.5;
delta = 1e-3;
Max_iter = 1000;

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
ttl = ['eta=' num2str(eta) '_a=' num2str(a) '_gamma=' num2str(gamma) '_L=' num2str(L) '_N=' num2str(N) '.mat']%file nya harus .mat tipenya + variabel" diatas

X0 = [real(E0);imag(E0)];

Ei = [0.01 0.011];

[Norm,Ei,lamb,temp_eig,jacobian] = LLE_continuation_fourdif(X0,Ei,eta,gamma,a,N,D2,theta,delta,Max_iter)


%save(ttl)
toc