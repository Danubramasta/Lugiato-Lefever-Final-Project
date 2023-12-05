function F = LLE_pseudo_fourdif(Y,eta,gamma,a,N,D2,theta,delta,X0,Ei0)
X = Y(1:2*N);
Ei = Y(2*N+1);

F1 = LLE_fourdif(X,Ei,eta,gamma,a,N,D2);
F2 = theta*norm(X-X0)^2+(1-theta)*abs(Ei-Ei0)^2-delta;

F = [F1;F2];