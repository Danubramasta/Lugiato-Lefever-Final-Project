clc; clear all; close all;
tic
L = 10;

N = 200;
dx = L/N;

m = 2;
[x,D2] = fourdif(N,m);
D2 = D2*(2*pi/L)^m; % space scalling
x = [0:dx:(L-dx)]';
% D2 = laplacian_1D(N)/dx^2;

% E0 = ones(N,1);
% eta = 1; a = 1;
% gamma = 1;
% 
% ttl = ['fourdif_amplitude_eta=' num2str(eta) '_a=' num2str(a) '_gamma=' num2str(gamma) '_L=' num2str(L) '_N=' num2str(N) '.mat']
% Ei_all = 1.02:0.01:3;
% idx1 = 1;
% for Ei = Ei_all;
%     
%     X0 = [real(E0);imag(E0)];
%     % F = LLE(X0,Ei,eta,theta,a,N,k);
%     
%     options = optimset('Display','iter','algorithm','levenberg-marquardt');   % Option to display output
%     [X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei,eta,gamma,a,N,D2),X0,options);
%     E = X(1:N)+1i*X(N+1:2*N);
%     
%     J = jacobian_LLE_fourdif(X,eta,gamma,a,N,D2);
%     
%     % temp_eig2 = eig(jacobian);
%     [VV,DD] = eig(J);
%     temp_eig = diag(DD);
%     [lamb,id] = max(real(temp_eig));
%     % plot(real(temp_eig),imag(temp_eig),'o',real(temp_eig2),imag(temp_eig2),'*')
%     % xlim([-2 0])
%     
%     t0 = 0;
%     dt = 0.0001;
%     T = 200;
%     t = t0:dt:T;
%     Max_iter = length(t);
%     X = X+(1e-4)*VV(:,id);
%     % X = X+(1e-4)*rand(2*N,1);
%     
%     for idx = 1:Max_iter
%         k1 = dt*LLE_fourdif(X,Ei,eta,gamma,a,N,D2);
%         k2 = dt*LLE_fourdif(X+0.5*k1,Ei,eta,gamma,a,N,D2);
%         k3 = dt*LLE_fourdif(X+0.5*k2,Ei,eta,gamma,a,N,D2);
%         k4 = dt*LLE_fourdif(X+k3,Ei,eta,gamma,a,N,D2);
%             
%         X  = X+(k1+2*k2+2*k3+k4)/6;
%         E = X(1:N)+1i*X(N+1:2*N);
%         
%             if mod(idx,100000) == 0
% %                 figure(1)
% %                 plot(x,abs(E))
% %                 title(['t=' num2str(t(idx))])
% %                 getframe;
%         
%             end
%     end
%     Ei_max(idx1) = max(abs(E))';
%     Ei_min(idx1) = min(abs(E))';
%     
%     
%     idx1 = idx1+1;
% end
% 
% %save(ttl)
% toc