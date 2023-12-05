%% Test Runge kutta
clc; clear all; close all;
tic
L = 50;

N = 200;
dx = L/N;

m = 2;
[x,D2] = fourdif(N,m);
D2 = D2*(2*pi/L)^m; % space scalling
x = [0:dx:(L-dx)]';


E0 = ones(N,1);
Ei = 1.5;
EE = [];
eta = 1; a = 1;
gamma = 1;
X0 = [real(E0);imag(E0)];
% F = LLE(X0,Ei,eta,theta,a,N,k);

options = optimset('Display','iter','algorithm','levenberg-marquardt');   % Option to display output
[X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei,eta,gamma,a,N,D2),X0,options);
E = X(1:N)+1i*X(N+1:2*N);

J = jacobian_LLE_fourdif(X,eta,gamma,a,N,D2);

[VV,DD] = eig(J); %ini vektor Eigen & matriks diagonal (isinya eigen)
temp_eig = diag(DD); %diambil diagonal matriks nya, makanya jadi vektor kolom yg isinya nilai eigen
[lamb,id] = max(real(temp_eig)); %didapat nilai eigen max di idx 368

t0 = 0;
dt = 0.001;
T = 200;
t = t0:dt:T;
Max_iter = length(t);
X = X+(1e-4)*VV(:,id);
A = LLE_fourdif(X,Ei,eta,gamma,a,N,D2);
B = LLE_fourdif(X+0.5*A,Ei,eta,gamma,a,N,D2);
% for idx = 1
%     k1 = dt*LLE_fourdif(X,Ei,eta,gamma,a,N,D2);
%     k2 = dt*LLE_fourdif(X+0.5*k1,Ei,eta,gamma,a,N,D2);
%     k3 = dt*LLE_fourdif(X+0.5*k2,Ei,eta,gamma,a,N,D2);
%     k4 = dt*LLE_fourdif(X+k3,Ei,eta,gamma,a,N,D2);
%     
%     X  = X+(k1+2*k2+2*k3+k4)/6;
%     E = X(1:N)+1i*X(N+1:2*N);
%     
%     EA = abs(E);
%     EN = min(abs(E));
%     EX = max(abs(E));
%     EE = [min(abs(E)) max(abs(E))];
% end
toc