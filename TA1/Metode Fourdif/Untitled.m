clc;clear all; close all;

% %% Tes Numerical Continuation
% theta = 0.5;
% delta = 1e-3;
% Max_iter = 1;
% 
% L = 50; %gaperlu phi karena phi untuk fungsi yang periodik sedangkan sekarang fungsi nya datar
% N = 200;
% dx = L/N;
% 
% %kebalik woy!!!
% m = 2;
% [x, D2] = fourdif(N,m); %dengan N banyaknya titik dan m turunan ke berapa
% D2 = D2*(2*pi/L)^m; %dinamakan space scalling
% x = [0:dx:(L-dx)]';
% 
% E0 = ones(N,1);
% 
% eta = 1; a = 1;
% gamma = 1;
% X0 = [real(E0);imag(E0)];
% 
% Ei = [0.01 0.011];
% options = optimset('Display','off','algorithm','levenberg-marquardt');   % Option to display output
% for idx = 1:3
%     if idx == 1
%         [X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei(idx),eta,gamma,a,N,D2),X0,options);
%         X1 = X;
%     elseif idx == 2
%         [X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei(idx),eta,gamma,a,N,D2),X1,options);
%         X2 = X;
%     else
%         X3 = 2*X2-X1;
%         Ei3 = 2*Ei(idx-1)-Ei(idx-2);
%         Y0 = [X3;Ei3];
%         [Y,fval,exitflag,output,jacobian] = fsolve(@(Y)LLE_pseudo_fourdif(Y,eta,gamma,a,N,D2,theta,delta,X2,Ei(idx-1)),Y0,options);
%         X = Y(1:2*N);
%         Ei(idx) = Y(2*N+1);
%         X1 = X2; X2 = X;
%     end
% end
% [X1 X2 X3]
% %% Newton Rapson
% eta = 1; a = 1;
% gamma = 1;
% N = 200;
% E0 = ones(N,1);
% X0 = [real(E0);imag(E0)];
% 
% options = optimset('Display','off','algorithm','levenberg-marquardt');   % Option to display output
% 
% [X,fval,exitflag] = fsolve(@(X)LLE_fourdif(X,Ei,eta,gamma,a,N,D2),X0,options);
%%
L = 10;

N = 200;
dx = L/N;

m = 2;
P =(2*pi/L)^m;
[x,D2] = fourdif(N,m);
%D2 = D2*(2*pi/L)^m; % space scalling
x = [0:dx:(L-dx)]';

N = 200;
h=2*pi/N;
b = -pi^2/3/h^2+1/12;
a= -pi^2/3/h^2-1/6;

n1=floor((N-1)/2); % pembulatan ke paling kecil
n2=ceil((N-1)/2); %pembulatan ke paling besar
kk=(1:N-1)';
topc=csc((1:n2)'*h/2).^2;
top1=csc((1*h/2)).^2;
top2=csc((2*h/2)).^2;
top3=csc((3*h/2)).^2;
c= -0.5*((-1).^kk);
d=-0.5*((-1).^kk).*[topc; flipud(topc(1:n1))];
m32 = -0.5*top1;
m42 = -(-1)^2*0.5*top2;
col1=[-pi^2/3/h^2-1/6; -0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];

% %m==2,                             % compute first column
% if rem(N,2)==0                         % of 2nd derivative matrix
%     topc=csc((1:n2)'*h/2).^2;
%     col1=[-pi^2/3/h^2-1/6; -0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
% else
%     topc=csc((1:n2)'*h/2).*cot((1:n2)'*h/2);
%     col1=[-pi^2/3/h^2+1/12; -0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];
% end;
% row1=col1;
% DM=toeplitz(col1,row1);