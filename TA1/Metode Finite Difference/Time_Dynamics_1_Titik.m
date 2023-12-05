clc; clear all; close all;
tic
L = 10;

N = 200;
dx = L/N;

x = [0:dx:(L-dx)]';
D2 = laplacian_1D(N)/dx^2;

E0 = ones(N,1);
eta = 1; a = 1;
gamma = 1;
Ei = 1.5;
X0 = [real(E0);imag(E0)];

options = optimset('Display','iter','algorithm','levenberg-marquardt');   % Option to display output
[X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei,eta,gamma,a,N,D2),X0,options);
E = X(1:N)+1i*X(N+1:2*N);

J = jacobian_LLE_fourdif(X,eta,gamma,a,N,D2);

[VV,DD] = eig(J);
temp_eig = diag(DD);
[lamb,id] = max(real(temp_eig));

t0 = 0;
dt = 0.0001;
T = 200;
t = t0:dt:T;
Max_iter = length(t);
X = X+(1e-4)*VV(:,id);

idx1 = 1;
filename = ['Time Dynamics untuk Ei =' num2str(Ei)]
for idx = 1:Max_iter
    k1 = dt*LLE_fourdif(X,Ei,eta,gamma,a,N,D2);
    k2 = dt*LLE_fourdif(X+0.5*k1,Ei,eta,gamma,a,N,D2);
    k3 = dt*LLE_fourdif(X+0.5*k2,Ei,eta,gamma,a,N,D2);
    k4 = dt*LLE_fourdif(X+k3,Ei,eta,gamma,a,N,D2);
    
    X  = X+(k1+2*k2+2*k3+k4)/6;
    E = X(1:N)+1i*X(N+1:2*N);
    
    if mod(idx,10000) == 0
        figure(1)
        plot(x,abs(E))
        title(['Time Dynamics untuk Ei =' num2str(Ei) ', t=' num2str(t(idx))])
        getframe;
        F(idx1) = getframe(gcf);
        idx1 = idx1 + 1;
    end
    Ei_max = max(abs(E));
    Ei_min = min(abs(E));
    EE = [Ei_max Ei_min];
end
EE
video = VideoWriter([filename '.mp4'],'MPEG-4');
video.FrameRate = 10;%4;
video.Quality = 100;
open(video)
writeVideo(video,F)
close(video)
clear video
save([filename '.mat'])
toc