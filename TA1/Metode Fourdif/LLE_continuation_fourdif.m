function [Norm,Ei,lamb,X,jacobian] = LLE_continuation_fourdif(X0,Ei,eta,gamma,a,N,D2,theta,delta,Max_iter)

options = optimset('Display','off','algorithm','levenberg-marquardt');   % Option to display output
for idx = 1:Max_iter
    if idx == 1
        [X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei(idx),eta,gamma,a,N,D2),X0,options);
        X1 = X;
    elseif idx == 2
        [X,fval,exitflag,output,jacobian] = fsolve(@(X)LLE_fourdif(X,Ei(idx),eta,gamma,a,N,D2),X1,options);
        X2 = X;
    else
        X3 = 2*X2-X1;
        Ei3 = 2*Ei(idx-1)-Ei(idx-2);
        Y0 = [X3;Ei3];
        [Y,fval,exitflag,output,jacobian] = fsolve(@(Y)LLE_pseudo_fourdif(Y,eta,gamma,a,N,D2,theta,delta,X2,Ei(idx-1)),Y0,options);
        X = Y(1:2*N);
        Ei(idx) = Y(2*N+1);
        X1 = X2; X2 = X;
    end
    E = X(1:N)+1i*X(N+1:2*N);
    [Ei(idx) abs(E(1))]
    Norm(idx) = abs(E(1));
    J = jacobian_LLE_fourdif(X,eta,gamma,a,N,D2);
    temp_eig = eig(J(1:2*N,1:2*N)); %nilai eigen
    lamb(idx) = max(real(temp_eig));

end