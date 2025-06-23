% This code is associated with the paper 
% "Unconditionally stable space–time isogeometric method for the linear Schrödinger equation"
%
% by M. Ferrari and S. Gómez

% This code reproduces Figure 1

clc
clear
close all

T = 1;
size_rho = 1000;
rho_vec = linspace(1,100,size_rho);
n = 1000;

cond_K = zeros(size_rho,6);

for p = 1:6

    cont_rho = 0;
    Nt = n-p+1;

    knots_t = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';
    t = ((0:Nt)/Nt)*T;
    t = t(2:end);
    h = t(1);

    B = mat_splines(Nt,p,1,1,knots_t,0,T);
    C = mat_splines(Nt,p,0,1,knots_t,0,T);
    B = B*h;
    B = B(1:end-1,2:end);
    C = C(1:end-1,2:end);
    
    for rho = rho_vec
        
        cont_rho = cont_rho+1;
        K = 1i*B-rho*C;

        cond_K(cont_rho,p) = cond(K);

    end
end

figure(1)
plot(rho_vec,cond_K(:,1))
hold on
plot(rho_vec,cond_K(:,2))
plot(rho_vec,cond_K(:,3))
legend('p=1','p=2','p=3')

figure(2)
plot(rho_vec,cond_K(:,4))
hold on
plot(rho_vec,cond_K(:,5))
plot(rho_vec,cond_K(:,6))
legend('p=4','p=5','p=6')